"""Query the NASA Exoplanet Archive and rank tonight's observable transits."""
import csv
import io
from dataclasses import dataclass, field
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import requests
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

from seestar_core import calculate_max_exposure, snr_full
from seestar_site import DarkWindow, local_iso, min_altitude_deg

ARCHIVE_URL = "https://exoplanetarchive.ipac.caltech.edu/TAP/sync"
EARTH_RADIUS_IN_SUN = 0.009157683   # R_earth / R_sun, for deriving Rp/Rs where not given directly


def sfloat(v, default=None):
    """Safe string-to-float; returns default for empty / NaN values."""
    try:
        x = float(v)
        return x if np.isfinite(x) else default
    except (TypeError, ValueError):
        return default


def fetch_transiting_planets(bright_limit: float, faint_limit: float, min_depth_pct: float,
                              timeout: int = 45) -> List[dict]:
    """Query the NASA Exoplanet Archive (TAP/ADQL) for transiting planets in the given
    magnitude / depth range. `min_depth_pct` is in percent (e.g. 1.0 = 1%)."""
    adql = (
        "SELECT pl_name,hostname,ra,dec,sy_vmag,"
        "pl_orbper,pl_tranmid,pl_trandur,pl_trandep,pl_ratror,pl_imppar "
        "FROM pscomppars "
        "WHERE tran_flag=1 "
        "AND pl_orbper IS NOT NULL AND pl_tranmid IS NOT NULL "
        "AND pl_trandur IS NOT NULL AND pl_ratror IS NOT NULL "
        "AND pl_trandep IS NOT NULL AND sy_vmag IS NOT NULL "
        f"AND sy_vmag>={bright_limit} AND sy_vmag<={faint_limit} "
        f"AND pl_trandep>={min_depth_pct}"
    )
    resp = requests.get(ARCHIVE_URL, params={"query": adql, "format": "csv"}, timeout=timeout)
    resp.raise_for_status()
    return list(csv.DictReader(io.StringIO(resp.text)))


def fetch_toi_candidates(bright_limit: float, faint_limit: float, min_depth_pct: float,
                          dispositions=("PC", "APC"), timeout: int = 45) -> List[dict]:
    """Query the NASA Exoplanet Archive's TESS Objects of Interest (TOI) table -
    the same catalog ExoFOP-TESS serves, mirrored on the archive's TAP endpoint.

    Unlike `pscomppars`, these are *unverified candidates*: some are eventual
    false positives. `dispositions` defaults to ('PC', 'APC') - planet
    candidate / ambiguous planet candidate - excluding confirmed planets
    (already covered by `fetch_transiting_planets`), false positives (FP),
    and false alarms (FA).

    Returns rows reshaped into the same dict schema as
    `fetch_transiting_planets`, with an added "source" key, so both can be
    passed to `find_tonight_transits` together.
    """
    disp_list = ",".join(f"'{d}'" for d in dispositions)
    adql = (
        "SELECT toi,tid,ra,dec,st_tmag,"
        "pl_orbper,pl_tranmid,pl_trandurh,pl_trandep,pl_rade,st_rad,tfopwg_disp "
        "FROM toi "
        "WHERE pl_orbper IS NOT NULL AND pl_tranmid IS NOT NULL "
        "AND pl_trandurh IS NOT NULL AND pl_trandep IS NOT NULL AND st_tmag IS NOT NULL "
        f"AND st_tmag>={bright_limit} AND st_tmag<={faint_limit} "
        f"AND pl_trandep>={min_depth_pct * 1e4} "   # pl_trandep is in ppm here (not % as in pscomppars)
        f"AND tfopwg_disp IN ({disp_list})"
    )
    resp = requests.get(ARCHIVE_URL, params={"query": adql, "format": "csv"}, timeout=timeout)
    resp.raise_for_status()
    rows = list(csv.DictReader(io.StringIO(resp.text)))

    planets = []
    for r in rows:
        dep_ppm = sfloat(r["pl_trandep"])
        if dep_ppm is None or dep_ppm <= 0:
            continue
        rade = sfloat(r["pl_rade"])
        srad = sfloat(r["st_rad"])
        if rade is not None and srad is not None and srad > 0:
            k = rade * EARTH_RADIUS_IN_SUN / srad          # Rp/Rs from radii
        else:
            k = np.sqrt(dep_ppm / 1e6)                      # fall back to depth = K^2
        if k <= 0:
            continue
        planets.append({
            "pl_name": f"TOI-{r['toi']}",
            "hostname": f"TIC {r['tid']}",
            "ra": r["ra"], "dec": r["dec"],
            "sy_vmag": r["st_tmag"],                         # TESS mag used as a V-mag proxy
            "pl_orbper": r["pl_orbper"],
            "pl_tranmid": r["pl_tranmid"],
            "pl_trandur": r["pl_trandurh"],
            "pl_trandep": str(dep_ppm / 1e4),                # ppm -> % to match pscomppars units
            "pl_ratror": str(k),
            "pl_imppar": "",                                 # not published for TOIs; defaults to 0 (central)
            "source": "TOI candidate",
        })
    return planets


def transits_tonight(t0_jd, per_d, dur_h, win0_jd, win1_jd):
    """JD mid-times of transits whose +/- half-duration overlaps [win0, win1]."""
    hd = dur_h / 48.0
    n0 = int(np.floor((win0_jd - hd - t0_jd) / per_d))
    n1 = int(np.ceil((win1_jd + hd - t0_jd) / per_d))
    return [t0_jd + n * per_d for n in range(n0, n1 + 1)
            if t0_jd + n * per_d + hd >= win0_jd
            and t0_jd + n * per_d - hd <= win1_jd]


@dataclass
class TransitCandidate:
    name: str
    host: str
    vmag: float
    coord: SkyCoord
    tmid_jd: float
    t1_jd: float
    t4_jd: float
    dur: float
    dep_pct: float
    k: float
    bimp: float
    alt: float
    frac: float
    exp: float
    snr5: float
    score: float
    source: str = "Confirmed"


def find_tonight_transits(planets: List[dict], dark_window: DarkWindow,
                           min_alt: float, sky_mag_obs: float, n_top: int = 5) -> List[TransitCandidate]:
    """Filter/rank tonight's observable transits from a planet list (see
    `fetch_transiting_planets`) against a computed dark window."""
    obs_loc = dark_window.obs_loc
    ds_jd, de_jd = dark_window.ds_jd, dark_window.de_jd

    candidates = []
    for p in planets:
        t0 = sfloat(p["pl_tranmid"]); per = sfloat(p["pl_orbper"])
        dur = sfloat(p["pl_trandur"]); dep = sfloat(p["pl_trandep"])
        vmag = sfloat(p["sy_vmag"]); ra = sfloat(p["ra"])
        dec = sfloat(p["dec"]); k = sfloat(p["pl_ratror"])
        bimp = sfloat(p["pl_imppar"], 0.0)

        if None in (t0, per, dur, dep, vmag, ra, dec, k) or k <= 0:
            continue
        bimp = np.clip(abs(bimp), 0.0, 1.0 + k - 1e-4)

        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")

        for mid_jd in transits_tonight(t0, per, dur, ds_jd, de_jd):
            t1j = mid_jd - dur / 48.0
            t4j = mid_jd + dur / 48.0
            ot0 = max(t1j, ds_jd)
            ot1 = min(t4j, de_jd)
            if ot1 <= ot0:
                continue
            alt = min_altitude_deg(coord, ot0, ot1, obs_loc)
            if alt < min_alt:
                continue

            frac = (ot1 - ot0) / (t4j - t1j)
            exp = max(1.0, min(calculate_max_exposure(vmag) * 0.8, 60.0))
            snr5 = snr_full(vmag, exp, sky_mag_arcsec2=sky_mag_obs) * np.sqrt(300.0 / exp)
            dep_pct = dep
            score = dep_pct * min(snr5, 1000) / 100 * frac * np.sin(np.radians(alt))

            candidates.append(TransitCandidate(
                name=p["pl_name"], host=p["hostname"], vmag=vmag, coord=coord,
                tmid_jd=mid_jd, t1_jd=t1j, t4_jd=t4j,
                dur=dur, dep_pct=dep_pct, k=k, bimp=bimp,
                alt=alt, frac=frac, exp=exp, snr5=snr5, score=score,
                source=p.get("source", "Confirmed"),
            ))

    seen = {}
    for c in candidates:
        if c.name not in seen or c.score > seen[c.name].score:
            seen[c.name] = c
    top = sorted(seen.values(), key=lambda x: -x.score)[:n_top]
    return top


def build_transit_figures(top: List[TransitCandidate], dark_window: DarkWindow,
                           sky_mag_obs: float, bin_min: float = 5.0, seed: int = 42) -> Optional[plt.Figure]:
    """Simulate + plot predicted light curves for the ranked candidates."""
    if not top:
        return None

    ds_jd, de_jd = dark_window.ds_jd, dark_window.de_jd

    ng = 160
    xg = np.linspace(-1.0, 1.0, ng)
    Xg, Yg = np.meshgrid(xg, xg)
    Rg = np.sqrt(Xg ** 2 + Yg ** 2)
    mu = np.where(Rg <= 1, np.sqrt(np.clip(1 - Rg ** 2, 0, None)), 0.0)
    u1d, u2d = 0.40, 0.27
    I_d = np.where(Rg <= 1, 1 - u1d * (1 - mu) - u2d * (1 - mu) ** 2, 0.0)
    I_dt = I_d.sum()

    rng = np.random.default_rng(seed)
    n_tp = len(top)
    fig, axes = plt.subplots(n_tp, 1, figsize=(14, 3.6 * n_tp), squeeze=False)
    axes = axes.ravel()
    fig.suptitle(
        f"Predicted Transits Tonight - {dark_window.site.name}\n"
        f"Dark window: {dark_window.dark_start.utc.iso[:16]} -> {dark_window.dark_end.utc.iso[:16]} UTC  "
        f"({local_iso(dark_window.dark_start, dark_window.site.tz)} -> "
        f"{local_iso(dark_window.dark_end, dark_window.site.tz)} {dark_window.site.tz})",
        fontsize=11, fontweight="bold"
    )

    for ax, c in zip(axes, top):
        k, bimp, dur, vmag, exp = c.k, c.bimp, c.dur, c.vmag, c.exp

        dn = np.sqrt(max(1e-9, (1 + k) ** 2 - bimp ** 2))
        dn23 = max(0.0, (1 - k) ** 2 - bimp ** 2)
        t23h = dur * np.sqrt(dn23) / dn

        obs_h = dur + 1.0
        n_f = int(obs_h * 3600 / exp) + 1
        t_h = np.linspace(-obs_h / 2, obs_h / 2, n_f)
        t_m = t_h * 60

        z_arr = np.sqrt(bimp ** 2 + (2.0 * t_h / dur * dn) ** 2)
        flux_m = np.ones(n_f)
        batch = 50
        for i0 in range(0, n_f, batch):
            zb = z_arr[i0:i0 + batch]
            if np.all(zb >= 1.0 + k):
                continue
            d2 = (Xg[:, :, None] - zb[None, None, :]) ** 2 + Yg[:, :, None] ** 2
            blk = np.sum(I_d[:, :, None] * (d2 <= k ** 2), axis=(0, 1))
            flux_m[i0:i0 + batch] = 1.0 - blk / I_dt

        sigma = 1.0 / snr_full(vmag, exp, sky_mag_arcsec2=sky_mag_obs)
        flux_o = flux_m + rng.normal(0.0, sigma, n_f)

        b_edges = np.arange(t_m[0], t_m[-1] + bin_min, bin_min)
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
                lbl = "Twilight / day" if first else "_"
                ax.axvspan(lo_m, hi_m, color="#eeeeee", zorder=0, label=lbl)
                first = False

        ax.errorbar(t_m, flux_o, yerr=sigma, fmt=".", ms=1.5,
                     color="0.70", elinewidth=0.2, alpha=0.22, zorder=1)
        ax.errorbar(b_ctrs, bfl, yerr=ber, fmt="o", ms=4.5, color="steelblue",
                     capsize=2.5, elinewidth=1.1, zorder=3, label=f"{bin_min:.0f}-min bins")
        ax.plot(t_m, flux_m, "r-", lw=2.0, zorder=4, label="Model")

        for tc, lab in [(-dur / 2 * 60, "T1"), (dur / 2 * 60, "T4")]:
            ax.axvline(tc, color="tomato", lw=0.8, ls="--", alpha=0.6)
            ax.text(tc, 1 - k ** 2 * 1.6, lab, ha="center", va="top", fontsize=7.5, color="tomato")
        if t23h > 0:
            for tc, lab in [(-t23h / 2 * 60, "T2"), (t23h / 2 * 60, "T3")]:
                ax.axvline(tc, color="tomato", lw=0.8, ls="--", alpha=0.6)
                ax.text(tc, 1 - k ** 2 * 1.6, lab, ha="center", va="top", fontsize=7.5, color="tomato")

        mid_str = Time(mid_jd, format="jd", scale="tdb").utc.iso[:16]
        src_tag = "" if c.source == "Confirmed" else f"  [{c.source}]"
        ax.set_title(
            f"{c.name}{src_tag}   V={vmag:.1f}   depth={c.dep_pct:.2f}%   "
            f"dur={dur:.1f} h   mid-transit {mid_str} UTC   "
            f"alt>={c.alt:.0f} deg   coverage {c.frac * 100:.0f}%   "
            f"SNR/5min~={c.snr5:.0f}   t_exp={exp:.0f} s",
            fontsize=8.5
        )
        ax.set_ylabel("Rel. flux", fontsize=9)
        ax.set_xlim(t_m[0], t_m[-1])
        ax.grid(alpha=0.2)
        if ax is axes[0]:
            ax.legend(fontsize=8, loc="lower right", ncol=3)

    axes[-1].set_xlabel("Time from mid-transit (min)", fontsize=10)
    plt.tight_layout()
    return fig
