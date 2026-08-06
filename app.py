"""
Seestar Exoplanet Transit Observation Planner - Streamlit app.

Four tools, one sidebar:
  1. SNR Explorer            - SNR vs exposure / magnitude / sky brightness
  2. Light Curve Simulator   - simulated transit for a single target
  3. Tonight's Best Transits - NASA Exoplanet Archive, ranked for your site
  4. Tonight's Best Eclipsing Binaries - AAVSO VSX + curated compact objects

Run with:  streamlit run app.py
"""
import matplotlib
matplotlib.use("Agg")

import streamlit as st

from seestar_eb import build_eb_figures, find_tonight_eclipses, query_vsx, vsx_cone_for_site
from seestar_lightcurve import simulate_transit
from seestar_site import Site, compute_dark_window, is_valid_tz, local_iso
from seestar_snr import build_snr_overview_figure
from seestar_transits import (
    build_transit_figures,
    fetch_toi_candidates,
    fetch_transiting_planets,
    find_tonight_transits,
)

st.set_page_config(page_title="Seestar Exoplanet Planner", layout="wide")

DEFAULT_SITE = Site(name="Frankfurt am Main, Germany", lat=50.0853, lon=8.5822, elev_m=100,
                     tz="Europe/Berlin")

if "site" not in st.session_state:
    st.session_state.site = DEFAULT_SITE


# =============================================================================
# Cached network calls
# =============================================================================
@st.cache_data(ttl=3600, show_spinner="Querying NASA Exoplanet Archive...")
def cached_fetch_transiting_planets(bright_limit, faint_limit, min_depth_pct):
    return fetch_transiting_planets(bright_limit, faint_limit, min_depth_pct)


@st.cache_data(ttl=3600, show_spinner="Querying TESS Objects of Interest (ExoFOP/TOI)...")
def cached_fetch_toi_candidates(bright_limit, faint_limit, min_depth_pct):
    return fetch_toi_candidates(bright_limit, faint_limit, min_depth_pct)


@st.cache_data(ttl=3600, show_spinner="Querying AAVSO VSX (this can take 1-2 min in dense fields)...")
def cached_query_vsx(ra_center, dec_center, radius_deg, bright_limit, faint_limit):
    return query_vsx(ra_center, dec_center, radius_deg, bright_limit, faint_limit)


# =============================================================================
# Site input widget, shared by pages 3 & 4
# =============================================================================
def site_form(key_prefix: str) -> Site:
    site = st.session_state.site
    with st.expander("Observer site", expanded=True):
        c1, c2, c3, c4, c5 = st.columns([2, 1, 1, 1, 1.4])
        name = c1.text_input("Site name", value=site.name, key=f"{key_prefix}_name")
        lat = c2.number_input("Latitude (deg N)", value=site.lat, min_value=-90.0, max_value=90.0,
                               format="%.4f", key=f"{key_prefix}_lat")
        lon = c3.number_input("Longitude (deg E)", value=site.lon, min_value=-180.0, max_value=180.0,
                               format="%.4f", key=f"{key_prefix}_lon")
        elev = c4.number_input("Elevation (m)", value=float(site.elev_m), min_value=-500.0,
                                max_value=9000.0, key=f"{key_prefix}_elev")
        tz = c5.text_input("Timezone (IANA)", value=site.tz, key=f"{key_prefix}_tz",
                            help="e.g. Europe/Berlin, America/Los_Angeles, Asia/Tokyo. "
                                 "Used only to display local times alongside UTC.")
        if not is_valid_tz(tz):
            st.warning(f"Unrecognized timezone '{tz}' - falling back to UTC for local-time display.")
    new_site = Site(name=name, lat=lat, lon=lon, elev_m=elev, tz=tz)
    st.session_state.site = new_site
    return new_site


# =============================================================================
# Page 1: SNR Explorer
# =============================================================================
def page_snr_explorer():
    st.header("SNR Explorer")
    st.caption("SNR as a function of exposure time, star magnitude, and sky brightness "
               "(Seestar S50, 50 mm f/5, IMX462).")

    c1, c2, c3 = st.columns(3)
    ref_sky = c1.slider("Sky brightness for panels 1 & 2 (mag/arcsec2)", 16.0, 22.0, 20.5, 0.1)
    exp_fixed = c2.slider("Fixed exposure for panel 3 (s)", 1.0, 120.0, 10.0, 1.0)
    star_mags = c3.multiselect("Star magnitudes to overlay", list(range(6, 17)),
                                default=[8, 10, 12, 13, 14, 15])
    exp_times = st.multiselect("Exposure times to overlay (panel 2, s)",
                                [1, 2, 5, 10, 20, 30, 60, 90, 120],
                                default=[1, 5, 10, 30, 60, 120])

    star_mags = sorted(star_mags) or [10]
    exp_times = sorted(exp_times) or [10]

    fig = build_snr_overview_figure(ref_sky=ref_sky, star_mags=star_mags, exp_times=exp_times,
                                     exp_fixed=exp_fixed, star_mags3=star_mags)
    st.pyplot(fig)


# =============================================================================
# Page 2: Light Curve Simulator
# =============================================================================
LIGHTCURVE_PRESETS = {
    "Custom": None,
    "HAT-P-32 b": dict(star_mag=11.44, k=0.1508, b_imp=0.147, t14_h=3.07, u1=0.30, u2=0.28),
    "HD 189733 b": dict(star_mag=7.65, k=0.1504, b_imp=0.671, t14_h=1.83, u1=0.30, u2=0.28),
}


def page_lightcurve_simulator():
    st.header("Light Curve Simulator")
    st.caption("Simulate a single-target transit light curve for given star/planet "
               "parameters and Seestar exposure settings.")

    preset_name = st.selectbox("Preset", list(LIGHTCURVE_PRESETS.keys()))
    preset = LIGHTCURVE_PRESETS[preset_name] or {}

    c1, c2, c3 = st.columns(3)
    with c1:
        star_mag = st.number_input("Host star magnitude (V)", value=preset.get("star_mag", 11.44),
                                    format="%.2f")
        k = st.number_input("Planet/star radius ratio (Rp/Rs)", value=preset.get("k", 0.1508),
                             min_value=0.001, max_value=0.9, format="%.4f")
        b_imp = st.number_input("Impact parameter", value=preset.get("b_imp", 0.147),
                                 min_value=0.0, max_value=1.5, format="%.3f")
        t14_h = st.number_input("Transit duration T14 (h)", value=preset.get("t14_h", 3.07),
                                 min_value=0.05, format="%.2f")
    with c2:
        u1 = st.number_input("Limb-darkening u1", value=preset.get("u1", 0.30), format="%.2f")
        u2 = st.number_input("Limb-darkening u2", value=preset.get("u2", 0.28), format="%.2f")
        exposure_s = st.number_input("Exposure time (s)", value=30.0, min_value=0.1, format="%.1f")
        sky_mag = st.number_input("Sky brightness (mag/arcsec2)", value=18.0, format="%.1f")
    with c3:
        obs_hours = st.number_input("Observation window (h)", value=5.0, min_value=0.5, format="%.1f")
        bin_min = st.number_input("Bin width (min)", value=5.0, min_value=0.1, format="%.1f")
        seed = st.number_input("Random seed", value=42, step=1)
        randomize = st.checkbox("Randomize seed on each run")

    if randomize:
        import numpy as _np
        seed = int(_np.random.default_rng().integers(0, 1_000_000))

    if k >= 1 - b_imp and b_imp < 1:
        st.warning("This is a grazing transit (b >= 1-K); the flat bottom (T2-T3) will vanish.")

    if st.button("Simulate", type="primary"):
        result = simulate_transit(star_mag=star_mag, k=k, b_imp=b_imp, t14_h=t14_h,
                                   u1=u1, u2=u2, exposure_s=exposure_s, sky_mag=sky_mag,
                                   obs_hours=obs_hours, bin_min=bin_min, seed=int(seed))
        st.pyplot(result.fig)

        m1, m2, m3, m4 = st.columns(4)
        m1.metric("Depth", f"{result.depth_pct:.3f}%", f"{result.depth_ppm:.0f} ppm")
        m2.metric("SNR / frame", f"{result.snr_per_frame:.1f}")
        m3.metric("Depth / frame sigma", f"{result.depth_over_sigma:.1f}x")
        m4.metric("Max safe exposure", f"{result.max_safe_exposure_s:.0f} s")
        st.caption(
            f"T14 = {result.t14_h:.2f} h, T23 = {result.t23_h:.2f} h, "
            f"ingress/egress {result.ingress_min:.1f} min each  |  "
            f"{result.frames_in_transit} frames -> {result.bins_in_transit} bins in transit, "
            f"binned sigma {result.binned_sigma_mmag:.2f} mmag/bin"
            + ("  |  GRAZING TRANSIT" if result.grazing else "")
        )


# =============================================================================
# Page 3: Tonight's Best Transits
# =============================================================================
def page_tonight_transits():
    st.header("Tonight's Best Observable Transits")
    st.caption("Queries the NASA Exoplanet Archive live and ranks transits observable "
               "from your site during tonight's dark window.")

    site = site_form("transits")

    c1, c2, c3 = st.columns(3)
    min_alt = c1.slider("Min. altitude throughout transit (deg)", 0.0, 80.0, 35.0, 1.0)
    bright_limit = c1.number_input("Bright limit (V mag, may saturate below)", value=5.5, format="%.1f")
    faint_limit = c2.number_input("Faint limit (V mag, poor SNR above)", value=10.5, format="%.1f")
    min_depth = c2.number_input("Min. transit depth (%)", value=1.5, min_value=0.0, format="%.2f")
    sky_mag_obs = c3.number_input("Sky brightness at your site (mag/arcsec2)", value=18.0, format="%.1f")
    n_top = c3.slider("Number of candidates to show", 1, 15, 5)

    include_toi = st.checkbox(
        "Include unverified TOI candidates (TESS Objects of Interest)",
        value=False,
        help="Adds TESS Objects of Interest not yet confirmed as planets (same catalog "
             "ExoFOP-TESS serves). Excludes known false positives/false alarms, but depths "
             "and durations are less reliable than for confirmed planets - some will turn "
             "out to be false positives.",
    )

    if st.button("Search tonight's transits", type="primary"):
        dark_window = compute_dark_window(site)
        if dark_window.placeholder:
            st.warning("No true dark window found in the next 36 h at this site - "
                       "using an 8-hour placeholder window.")
        st.write(
            f"**Site:** {site.name} ({site.lat:+.3f}, {site.lon:+.3f})  \n"
            f"**Dark window:** {dark_window.dark_start.utc.iso[:16]} - "
            f"{dark_window.dark_end.utc.iso[:16]} UTC "
            f"({(dark_window.de_jd - dark_window.ds_jd) * 24:.1f} h)  \n"
            f"**Local time ({site.tz}):** {local_iso(dark_window.dark_start, site.tz)} - "
            f"{local_iso(dark_window.dark_end, site.tz)}"
        )

        planets = cached_fetch_transiting_planets(bright_limit, faint_limit, min_depth)
        st.write(f"{len(planets)} confirmed transiting planets in magnitude/depth range retrieved from archive.")

        if include_toi:
            toi_planets = cached_fetch_toi_candidates(bright_limit, faint_limit, min_depth)
            st.write(f"{len(toi_planets)} unverified TOI candidates in range retrieved from ExoFOP/TOI table.")
            planets = planets + toi_planets

        top = find_tonight_transits(planets, dark_window, min_alt, sky_mag_obs, n_top)

        if not top:
            st.info("No transits found. Try relaxing min. altitude, faint limit, or min. depth.")
            return

        st.subheader(f"Top {len(top)} candidates")
        st.dataframe([{
            "Planet": c.name, "Host": c.host, "Source": c.source, "V": round(c.vmag, 1),
            "Depth %": round(c.dep_pct, 2), "Duration (h)": round(c.dur, 1),
            "Min alt (deg)": round(c.alt, 0), "Coverage %": round(c.frac * 100, 0),
            "SNR/5min": round(c.snr5, 0),
            "Mid-transit (UTC)": c.tmid_jd,
        } for c in top], hide_index=True)

        fig = build_transit_figures(top, dark_window, sky_mag_obs)
        if fig is not None:
            st.pyplot(fig)


# =============================================================================
# Page 4: Tonight's Best Eclipsing Binaries
# =============================================================================
def page_tonight_eb():
    st.header("Tonight's Best Eclipsing Binaries & Compact-Object Systems")
    st.caption("Queries AAVSO VSX live (cone search around zenith) plus a curated list of "
               "compact-object eclipsing systems (WD/NS/BH), ranked for your site tonight. "
               "The VSX query can take 1-2 minutes in dense fields.")

    site = site_form("eb")

    c1, c2, c3 = st.columns(3)
    min_depth_eb = c1.number_input("Min. primary eclipse depth (mag)", value=0.01, min_value=0.0, format="%.3f")
    max_period_eb = c1.number_input("Max. period (days)", value=365.0, min_value=0.01, format="%.1f")
    min_dur_eb_min = c2.number_input("Min. eclipse duration (min)", value=5.0, min_value=0.1, format="%.1f")
    max_dur_eb = c2.number_input("Max. eclipse duration (h)", value=2.0, min_value=0.1, format="%.2f")
    min_alt_eb = c3.slider("Min. altitude throughout eclipse (deg)", 0.0, 80.0, 50.0, 1.0,
                            help="Kept high by default: AAVSO VSX reliably times out above "
                                 "~45-50 deg cones on dense fields.")
    bright_limit_eb = c1.number_input("Bright limit (V mag)", value=5.0, format="%.1f")
    faint_limit_eb = c2.number_input("Faint limit (V mag)", value=13.0, format="%.1f")
    sky_mag_obs_eb = c3.number_input("Sky brightness at your site (mag/arcsec2)", value=18.0, format="%.1f")
    n_top_eb = st.slider("Number of candidates to show", 1, 20, 10)
    include_compact = st.checkbox("Include curated compact-object systems (WD/NS/BH)", value=True)

    if st.button("Search tonight's eclipsing binaries", type="primary"):
        dark_window = compute_dark_window(site)
        if dark_window.placeholder:
            st.warning("No true dark window found in the next 36 h at this site - "
                       "using an 8-hour placeholder window.")
        st.write(
            f"**Site:** {site.name} ({site.lat:+.3f}, {site.lon:+.3f})  \n"
            f"**Dark window:** {dark_window.dark_start.utc.iso[:16]} - "
            f"{dark_window.dark_end.utc.iso[:16]} UTC "
            f"({(dark_window.de_jd - dark_window.ds_jd) * 24:.1f} h)  \n"
            f"**Local time ({site.tz}):** {local_iso(dark_window.dark_start, site.tz)} - "
            f"{local_iso(dark_window.dark_end, site.tz)}"
        )

        ra_center, dec_center, radius_deg = vsx_cone_for_site(dark_window, min_alt_eb)
        vsx_rows, used_radius = cached_query_vsx(ra_center, dec_center, radius_deg,
                                                  bright_limit_eb, faint_limit_eb)
        if used_radius is None:
            st.warning("AAVSO VSX unreachable even at the smallest retry radius - "
                       "falling back to the curated compact-object list only.")
        elif used_radius < radius_deg - 0.5:
            st.info(f"Had to shrink the search cone to r={used_radius:.0f} deg "
                    f"(target was {radius_deg:.0f} deg) - AAVSO couldn't finish the full "
                    f"query tonight, so coverage is reduced.")
        st.write(f"{len(vsx_rows)} eclipsing-type variables retrieved from VSX.")

        top = find_tonight_eclipses(
            vsx_rows, dark_window,
            bright_limit=bright_limit_eb, faint_limit=faint_limit_eb,
            min_depth_eb=min_depth_eb, max_period_eb=max_period_eb,
            min_dur_eb=min_dur_eb_min / 60.0, max_dur_eb=max_dur_eb,
            min_alt_eb=min_alt_eb, sky_mag_obs=sky_mag_obs_eb, n_top=n_top_eb,
            include_compact_systems=include_compact,
        )

        if not top:
            st.info("No eclipsing systems found. Try relaxing min. altitude, faint limit, or min. depth.")
            return

        st.subheader(f"Top {len(top)} candidates")
        st.dataframe([{
            "System": c.name, "Type": c.label, "V": round(c.vmag, 1),
            "Depth (mag)": round(c.depth_mag, 2), "Duration (h)": round(c.dur_h, 2),
            "Period (d)": round(c.per, 3),
            "Min alt (deg)": round(c.alt, 0), "Coverage %": round(c.frac * 100, 0),
            "SNR/5min": round(c.snr5, 0),
            "Mid-eclipse (UTC)": c.tmid_jd,
        } for c in top], hide_index=True)

        fig = build_eb_figures(top, dark_window, sky_mag_obs_eb)
        if fig is not None:
            st.pyplot(fig)


# =============================================================================
# Router
# =============================================================================
PAGES = {
    "SNR Explorer": page_snr_explorer,
    "Light Curve Simulator": page_lightcurve_simulator,
    "Tonight's Best Transits": page_tonight_transits,
    "Tonight's Best Eclipsing Binaries": page_tonight_eb,
}

st.sidebar.title("Seestar Exoplanet Planner")
choice = st.sidebar.radio("Tool", list(PAGES.keys()))
PAGES[choice]()
