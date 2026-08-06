"""Query AAVSO VSX + a curated compact-object list, and rank tonight's
observable eclipsing binaries / compact-object eclipses."""
import re
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from typing import List, Optional, Tuple

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from seestar_core import calculate_max_exposure, snr_full
from seestar_site import DarkWindow, local_iso, min_altitude_deg

VSX_URL = "https://www.aavso.org/vsx/index.php"

# name, RA(deg), Dec(deg), V_mag, period(d), epoch(JD),
# depth_primary(mag), dur_frac (eclipse duration / period), type_label
COMPACT_SYSTEMS = [
    # White dwarf + main-sequence eclipsing (WD+MS)
    ("V471 Tau",      52.6858,  17.2441,  9.4, 0.5212, 2448612.4638, 0.55, 0.05, "WD+MS"),
    ("HW Vir",       186.5917,  -7.5175, 10.8, 0.1167, 2446929.0977, 1.10, 0.07, "sdB+MS"),
    ("HS 0705+6700", 107.2792,  66.9028, 14.7, 0.0957, 2451822.7406, 0.55, 0.06, "sdB+MS"),
    ("IK Peg",       327.3671,  19.1487,  6.1, 21.722, 2450259.0960, 0.15, 0.008, "WD+MS"),
    ("NN Ser",       148.1885,  12.7321, 16.5, 0.1301, 2447344.6230, 0.92, 0.02, "WD+MS"),
    # Neutron star X-ray binary
    ("HZ Her",       254.4575,  35.3422, 13.5, 1.7001, 2440946.4260, 1.50, 0.14, "NS+MS"),
    ("V1341 Cyg",    326.1713,  38.3214, 14.7, 9.8444, 2440788.2000, 0.30, 0.05, "NS+MS"),
    # Black-hole X-ray binaries
    ("HDE 226868",   299.5904,  35.2016,  8.9, 5.5999, 2441163.5290, 0.08, 0.10, "BH+OB"),
    ("KV UMa",       169.5454,  48.0367, 18.8, 0.1699, 2451773.3800, 0.05, 0.10, "BH+MS"),
    ("V1033 Sco",    253.4863, -39.8458, 17.0, 2.6217, 2449838.8000, 0.20, 0.06, "BH+MS"),
]

_MAG_RE = re.compile(r"\(?(-?\d+\.?\d*)\)?\s*([A-Za-z]*)")


def parse_mag(s):
    """'9.94 V' -> (9.94, 'V');  '(0.54) CV' -> (0.54, 'CV')."""
    m = _MAG_RE.match((s or "").strip())
    if not m:
        return None, None
    try:
        return float(m.group(1)), m.group(2)
    except ValueError:
        return None, None


def is_eclipsing_type(vtype):
    tokens = re.split(r"[/|+:]", vtype or "")
    return any(tok.startswith(("EA", "EB", "EW")) for tok in tokens)


def vsx_cone_for_site(dark_window: DarkWindow, min_alt_eb: float) -> Tuple[float, float, float]:
    """RA/Dec/radius (deg) of the search cone centred on zenith at the
    midpoint of tonight's dark window, sized to reach `min_alt_eb`."""
    obs_loc = dark_window.obs_loc
    mid_jd = (dark_window.ds_jd + dark_window.de_jd) / 2.0
    mid_time = Time(mid_jd, format="jd", scale="tdb")
    ra_center = mid_time.sidereal_time("apparent", longitude=obs_loc.lon).deg % 360.0
    dec_center = obs_loc.lat.deg
    radius_deg = max(5.0, min(90.0, 90.0 - min_alt_eb))
    return ra_center, dec_center, radius_deg


def query_vsx(ra_center: float, dec_center: float, radius_deg: float,
              bright_limit: float, faint_limit: float,
              timeout: int = 150) -> Tuple[List[tuple], Optional[float]]:
    """Query AAVSO VSX for eclipsing-type variables (EA/EB/EW) in a cone.

    The VSX API only supports cone search (no server-side mag/type filtering,
    see notebook cell 10 notes) and always replies with XML regardless of the
    requested format. A dense field can time out on AAVSO's side above
    ~45-50 deg cones, so we shrink the radius and retry on failure.

    Returns (rows, used_radius_deg). Rows are
    (name, ra, dec, maxmag, period, epoch, amplitude, dur_frac, vtype) tuples,
    already filtered to eclipsing types + mag range + finite values.
    """
    root = None
    used_radius = None
    r_try = radius_deg
    while r_try >= 12.0:
        try:
            resp = requests.get(VSX_URL, params={"view": "api.list", "ra": ra_center,
                                                   "dec": dec_center, "radius": r_try},
                                 timeout=timeout)
            resp.raise_for_status()
            root = ET.fromstring(resp.content)
            used_radius = r_try
            break
        except (ET.ParseError, requests.RequestException):
            r_try *= 0.7

    rows = []
    if root is None:
        return rows, None

    for obj in root.findall("VSXObject"):
        vtype = obj.findtext("VariabilityType", "")
        if not is_eclipsing_type(vtype):
            continue
        try:
            ra = float(obj.findtext("RA2000", "nan"))
            dec = float(obj.findtext("Declination2000", "nan"))
            per = float(obj.findtext("Period", "nan"))
            ep = float(obj.findtext("Epoch", "nan"))
        except (TypeError, ValueError):
            continue
        maxmag, band1 = parse_mag(obj.findtext("MaxMag", ""))
        minmag, band2 = parse_mag(obj.findtext("MinMag", ""))
        if maxmag is None or minmag is None or band1 != band2:
            continue
        amp = minmag - maxmag
        if not all(np.isfinite(v) for v in (ra, dec, per, ep, maxmag, amp)):
            continue
        if per <= 0 or amp < 0:
            continue
        if maxmag < bright_limit or maxmag > faint_limit:
            continue
        try:
            dur_frac = float(obj.findtext("EclipseDuration", "")) / 100.0
            if not (0 < dur_frac < 1):
                dur_frac = 0.08
        except (TypeError, ValueError):
            dur_frac = 0.08
        rows.append((obj.findtext("Name", "?"), ra, dec, maxmag, per, ep, amp, dur_frac, vtype))

    return rows, used_radius


def eclipses_tonight(t0_jd, per_d, dur_frac, win0, win1):
    """JD mid-eclipse times overlapping the dark window [win0, win1]."""
    hd = per_d * dur_frac / 2.0
    n0 = int(np.floor((win0 - hd - t0_jd) / per_d))
    n1 = int(np.ceil((win1 + hd - t0_jd) / per_d))
    return [t0_jd + n * per_d for n in range(n0, n1 + 1)
            if t0_jd + n * per_d + hd >= win0
            and t0_jd + n * per_d - hd <= win1]


@dataclass
class EclipseCandidate:
    name: str
    label: str
    vmag: float
    coord: SkyCoord
    tmid_jd: float
    t1_jd: float
    t4_jd: float
    dur_h: float
    per: float
    depth_mag: float
    dep_frac: float
    alt: float
    frac: float
    exp: float
    snr5: float
    score: float


def find_tonight_eclipses(vsx_rows: List[tuple], dark_window: DarkWindow,
                           bright_limit: float, faint_limit: float,
                           min_depth_eb: float, max_period_eb: float,
                           min_dur_eb: float, max_dur_eb: float, min_alt_eb: float,
                           sky_mag_obs: float, n_top: int = 10,
                           include_compact_systems: bool = True) -> List[EclipseCandidate]:
    """Filter/rank tonight's observable eclipses from VSX rows + curated
    compact-object systems against a computed dark window."""
    obs_loc = dark_window.obs_loc
    ds_jd, de_jd = dark_window.ds_jd, dark_window.de_jd

    all_systems = list(vsx_rows)
    if include_compact_systems:
        all_systems += list(COMPACT_SYSTEMS)

    candidates = []
    for (name, ra, dec, vmag, per, ep_jd, depth_mag, dur_frac, label) in all_systems:
        if not all(np.isfinite(v) for v in (ra, dec, vmag, per, ep_jd, depth_mag)):
            continue
        if vmag < bright_limit or vmag > faint_limit:
            continue
        if per > max_period_eb or depth_mag < min_depth_eb:
            continue

        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
        dur_h = per * dur_frac * 24.0
        if not (min_dur_eb <= dur_h <= max_dur_eb):
            continue

        for mid_jd in eclipses_tonight(ep_jd, per, dur_frac, ds_jd, de_jd):
            t1j = mid_jd - per * dur_frac / 2.0
            t4j = mid_jd + per * dur_frac / 2.0
            ot0 = max(t1j, ds_jd)
            ot1 = min(t4j, de_jd)
            if ot1 <= ot0:
                continue
            alt = min_altitude_deg(coord, ot0, ot1, obs_loc)
            if alt < min_alt_eb:
                continue

            frac = (ot1 - ot0) / (t4j - t1j)
            exp = max(1.0, min(calculate_max_exposure(vmag) * 0.8, 60.0))
            snr5 = snr_full(vmag, exp, sky_mag_arcsec2=sky_mag_obs) * np.sqrt(300.0 / exp)
            dep_frac = 1.0 - 10 ** (-0.4 * depth_mag)
            score = dep_frac * min(snr5, 1000) / 100 * frac * np.sin(np.radians(alt))

            candidates.append(EclipseCandidate(
                name=name, label=label, vmag=vmag, coord=coord,
                tmid_jd=mid_jd, t1_jd=t1j, t4_jd=t4j,
                dur_h=dur_h, per=per, depth_mag=depth_mag, dep_frac=dep_frac,
                alt=alt, frac=frac, exp=exp, snr5=snr5, score=score,
            ))

    seen = {}
    for c in candidates:
        if c.name not in seen or c.score > seen[c.name].score:
            seen[c.name] = c
    top = sorted(seen.values(), key=lambda x: -x.score)[:n_top]
    return top


def trap_model(t_h_arr, half_dur_h, dep, ing_h):
    """Trapezoidal eclipse model centred at t=0."""
    flat_h = half_dur_h - ing_h
    flux = np.ones_like(t_h_arr)
    at = np.abs(t_h_arr)
    in_flat = at <= flat_h
    in_ingress = (at > flat_h) & (at <= flat_h + ing_h)
    flux[in_flat] = 1.0 - dep
    flux[in_ingress] = 1.0 - dep * (1.0 - (at[in_ingress] - flat_h) / ing_h)
    return flux


def build_eb_figures(top: List[EclipseCandidate], dark_window: DarkWindow,
                      sky_mag_obs: float, bin_min_eb: float = 0.5,
                      seed: int = 42) -> Optional[plt.Figure]:
    """Simulate + plot predicted (trapezoidal) eclipse curves for ranked candidates."""
    if not top:
        return None

    ds_jd, de_jd = dark_window.ds_jd, dark_window.de_jd
    rng_eb = np.random.default_rng(seed)
    n_eb = len(top)
    fig, axes = plt.subplots(n_eb, 1, figsize=(14, 3.4 * n_eb), squeeze=False)
    axes = axes.ravel()
    fig.suptitle(
        f"Predicted Eclipses Tonight - {dark_window.site.name}\n"
        f"Dark window: {dark_window.dark_start.utc.iso[:16]} -> {dark_window.dark_end.utc.iso[:16]} UTC  "
        f"({local_iso(dark_window.dark_start, dark_window.site.tz)} -> "
        f"{local_iso(dark_window.dark_end, dark_window.site.tz)} {dark_window.site.tz})",
        fontsize=11, fontweight="bold"
    )

    for ax, c in zip(axes, top):
        dur_h, exp, vmag, dep = c.dur_h, c.exp, c.vmag, c.dep_frac
        half_dur = dur_h / 2.0
        ing_h = half_dur * 0.15

        obs_h = dur_h + 1.0
        n_f = int(obs_h * 3600 / exp) + 1
        t_h = np.linspace(-obs_h / 2, obs_h / 2, n_f)
        t_m = t_h * 60.0

        flux_m = trap_model(t_h, half_dur, dep, ing_h)

        sigma = 1.0 / snr_full(vmag, exp, sky_mag_arcsec2=sky_mag_obs)
        flux_o = flux_m + rng_eb.normal(0.0, sigma, n_f)

        b_edges = np.arange(t_m[0], t_m[-1] + bin_min_eb, bin_min_eb)
        b_ctrs = (b_edges[:-1] + b_edges[1:]) / 2
        bfl, ber = [], []
        for lo, hi in zip(b_edges[:-1], b_edges[1:]):
            idx = (t_m >= lo) & (t_m < hi); nb = idx.sum()
            bfl.append(flux_o[idx].mean() if nb else np.nan)
            ber.append(sigma / np.sqrt(nb) if nb else np.nan)
        bfl = np.array(bfl); ber = np.array(ber)

        mid_jd = c.tmid_jd
        ds_m = (ds_jd - mid_jd) * 1440
        de_m = (de_jd - mid_jd) * 1440
        first = True
        for lo_m, hi_m in [(t_m[0], min(t_m[-1], ds_m)), (max(t_m[0], de_m), t_m[-1])]:
            if hi_m > lo_m:
                ax.axvspan(lo_m, hi_m, color="#eeeeee", zorder=0, label="Twilight / day" if first else "_")
                first = False

        ax.errorbar(t_m, flux_o, yerr=sigma, fmt=".", ms=1.5,
                     color="0.70", elinewidth=0.2, alpha=0.22, zorder=1)
        ax.errorbar(b_ctrs, bfl, yerr=ber, fmt="o", ms=4.5, color="steelblue",
                     capsize=2.5, elinewidth=1.1, zorder=3, label=f"{bin_min_eb:.0f}-min bins")
        ax.plot(t_m, flux_m, "r-", lw=2.0, zorder=4, label="Model (trapezoidal)")

        flat_m = (half_dur - ing_h) * 60
        tot_m = half_dur * 60
        for tc, lab in [(-tot_m, "I"), (-flat_m, "II"), (flat_m, "III"), (tot_m, "IV")]:
            ax.axvline(tc, color="tomato", lw=0.8, ls="--", alpha=0.6)
            ax.text(tc, 1 - dep * 1.6, lab, ha="center", va="top", fontsize=7.5, color="tomato")

        mid_str = Time(mid_jd, format="jd", scale="tdb").utc.iso[:16]
        ax.set_title(
            f"{c.name} [{c.label}]   V={vmag:.1f}   "
            f"depth={c.depth_mag:.2f} mag   dur={dur_h:.1f} h   "
            f"P={c.per:.4g} d   mid-eclipse {mid_str} UTC   "
            f"alt>={c.alt:.0f} deg   cov {c.frac * 100:.0f}%   "
            f"SNR/5min~={c.snr5:.0f}   t_exp={exp:.0f} s",
            fontsize=8.5
        )
        ax.set_ylabel("Rel. flux", fontsize=9)
        ax.set_xlim(t_m[0], t_m[-1])
        ax.grid(alpha=0.2)
        if ax is axes[0]:
            ax.legend(fontsize=8, loc="lower right", ncol=3)

    axes[-1].set_xlabel("Time from mid-eclipse (min)", fontsize=10)
    plt.tight_layout()
    return fig
