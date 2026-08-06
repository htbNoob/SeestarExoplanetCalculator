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

from seestar_core import SEESTAR_CAMERA, SEESTAR_OTA, Camera, Instrument, OTA
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
from seestar_weather import build_weather_figure, fetch_hourly_forecast, summarize_dark_window_weather

st.set_page_config(page_title="Seestar Exoplanet Planner", layout="wide")

MANUAL_OPTION = "Manual (custom location)"

PRESET_SITES = {
    "Frankfurt am Main, Germany": Site(
        name="Frankfurt am Main, Germany", lat=50.0853, lon=8.5822, elev_m=100,
        tz="Europe/Berlin", sky_mag=18.0),
    "Taunus Observatory (Kleiner Feldberg)": Site(
        name="Taunus Observatory (Kleiner Feldberg)", lat=50.2217, lon=8.4458, elev_m=826,
        tz="Europe/Berlin", sky_mag=20.8),    # Bortle 4
    "Teide Observatory, Tenerife": Site(
        name="Teide Observatory, Tenerife", lat=28.3000, lon=-16.5097, elev_m=2390,
        tz="Atlantic/Canary", sky_mag=21.59),  # Bortle 4
    "Madrona Peak Observatory": Site(
        name="Madrona Peak Observatory", lat=29.8442, lon=-99.3558, elev_m=379,
        tz="America/Chicago", sky_mag=21.89),  # Bortle 2
}
DEFAULT_PRESET = "Frankfurt am Main, Germany"
DEFAULT_SITE = PRESET_SITES[DEFAULT_PRESET]

if "site" not in st.session_state:
    st.session_state.site = DEFAULT_SITE

# =============================================================================
# Instrument presets (optical tube x camera, picked independently)
#
# Specs marked "est." below aren't in the source spec sheets and are
# estimated (from manufacturer literature or self-consistent from full
# well / bit depth) - worth confirming against your own calibration frames.
# =============================================================================
OTA_PRESETS = {
    "Seestar S50": SEESTAR_OTA,
    "Taunus Westkuppel (PlaneWave CDK12.5, f/8)": OTA(
        name="Taunus Westkuppel CDK12.5", aperture_mm=317.5, focal_length_mm=2531.0),
    "Taunus Ostkuppel 60cm (primary focus, f/3.3)": OTA(
        name="Taunus Ostkuppel 60cm (primary)", aperture_mm=600.0, focal_length_mm=2000.0),
    "Taunus Ostkuppel 60cm (secondary focus, f/10)": OTA(
        name="Taunus Ostkuppel 60cm (secondary)", aperture_mm=600.0, focal_length_mm=6000.0),
    "MPO61 RCOS 610 24\" (f/7.9)": OTA(
        name="MPO61 RCOS 610", aperture_mm=609.6, focal_length_mm=4814.0),
}

CAMERA_PRESETS = {
    "Seestar Sony IMX462": SEESTAR_CAMERA,
    "Atik 383L+ (KAF-8300 CCD)": Camera(
        name="Atik 383L+", pixel_size_um=5.4, full_well_e=25500, read_noise_e=8.5,
        dark_current_e_s=0.02, qe=0.5, gain_e_per_adu=0.39, bit_depth=16),
        # full_well/read_noise: midpoints of Atik's published 25000-26000e-/7-10e- ranges.
        # dark_current, gain: est. (cooled CCD typical / full_well / 2^16).
    "Moravian C3-61000 (Sony IMX455 CMOS)": Camera(
        name="Moravian C3-61000", pixel_size_um=3.76, full_well_e=50000, read_noise_e=2.0,
        dark_current_e_s=0.001, qe=0.85, gain_e_per_adu=0.76, bit_depth=16),
        # read_noise/qe: midpoints of Moravian's published 1-3e-/80-90% ranges.
        # dark_current, gain: est. (cooled BSI CMOS typical / full_well / 2^16).
    "QHY461PH (MPO61)": Camera(
        name="QHY461PH", pixel_size_um=3.76, full_well_e=75600, read_noise_e=3.5,
        dark_current_e_s=0.0005, qe=0.85, gain_e_per_adu=1.26, bit_depth=16),
        # full_well = linearitylimit(60000 ADU) x gain(1.26 e-/ADU) from MPO61's own config.
        # read_noise, dark_current, qe: est. (QHYCCD lit. quotes 1-3.7e- read noise by gain mode).
}
DEFAULT_OTA = "Seestar S50"
DEFAULT_CAMERA = "Seestar Sony IMX462"


def instrument_form(key_prefix: str) -> Instrument:
    with st.expander("Instrument", expanded=True):
        c1, c2 = st.columns(2)
        ota_name = c1.selectbox("Optical tube", list(OTA_PRESETS.keys()), key=f"{key_prefix}_ota")
        cam_name = c2.selectbox("Camera", list(CAMERA_PRESETS.keys()), key=f"{key_prefix}_cam")
        instrument = Instrument(ota=OTA_PRESETS[ota_name], camera=CAMERA_PRESETS[cam_name])
        st.caption(
            f"{instrument.ota.aperture_mm:.0f} mm aperture, f/{instrument.ota.f_ratio:.1f}  |  "
            f"{instrument.pixel_scale_arcsec:.2f}\"/px  |  QE {instrument.camera.qe * 100:.0f}%  |  "
            f"read noise {instrument.camera.read_noise_e:.1f} e-  |  "
            f"full well {instrument.camera.full_well_e:.0f} e-"
        )
    return instrument


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


@st.cache_data(ttl=1800, show_spinner="Fetching weather forecast (Open-Meteo)...")
def cached_fetch_hourly_forecast(lat, lon):
    return fetch_hourly_forecast(lat, lon)


def show_weather(dark_window):
    """Fetch + display the cloud cover / precipitation / temperature forecast
    for tonight's dark window. Shared by the transits and EB pages."""
    forecast = cached_fetch_hourly_forecast(dark_window.site.lat, dark_window.site.lon)
    summary = summarize_dark_window_weather(forecast, dark_window)
    st.subheader("Weather forecast")
    if summary is None:
        st.info("No forecast hours fall within tonight's dark window.")
        return
    m1, m2, m3, m4 = st.columns(4)
    m1.metric("Conditions", summary.rating)
    m2.metric("Mean cloud cover", f"{summary.mean_cloudcover:.0f}%", f"max {summary.max_cloudcover:.0f}%")
    m3.metric("Mean precip. probability", f"{summary.mean_precip_prob:.0f}%")
    m4.metric("Min. temp / max wind", f"{summary.min_temperature_c:.0f} C / {summary.max_windspeed_kmh:.0f} km/h")
    st.caption("Rating: Good = <20% mean cloud cover & <20% precip. probability during the dark window; "
               "Fair = <60% cloud & <40% precip.; otherwise Poor. Forecast: Open-Meteo, updated every 30 min.")
    st.pyplot(build_weather_figure(forecast, dark_window))


# =============================================================================
# Site input widget, shared by pages 3 & 4
# =============================================================================
def _apply_preset(key_prefix: str, preset_key: str):
    choice = st.session_state[preset_key]
    if choice == MANUAL_OPTION:
        return
    preset = PRESET_SITES[choice]
    st.session_state[f"{key_prefix}_name"] = preset.name
    st.session_state[f"{key_prefix}_lat"] = preset.lat
    st.session_state[f"{key_prefix}_lon"] = preset.lon
    st.session_state[f"{key_prefix}_elev"] = float(preset.elev_m)
    st.session_state[f"{key_prefix}_tz"] = preset.tz
    st.session_state[f"{key_prefix}_sky_mag"] = preset.sky_mag


def site_form(key_prefix: str) -> Site:
    site = st.session_state.site
    preset_key = f"{key_prefix}_preset"
    options = list(PRESET_SITES.keys()) + [MANUAL_OPTION]

    # Seed each widget's session-state value once; afterwards the key alone owns the
    # value (via direct edits or _apply_preset), so widgets below must not also pass
    # `value=` - doing both trips Streamlit's "value set via Session State API" warning.
    st.session_state.setdefault(f"{key_prefix}_name", site.name)
    st.session_state.setdefault(f"{key_prefix}_lat", site.lat)
    st.session_state.setdefault(f"{key_prefix}_lon", site.lon)
    st.session_state.setdefault(f"{key_prefix}_elev", float(site.elev_m))
    st.session_state.setdefault(f"{key_prefix}_tz", site.tz)

    with st.expander("Observer site", expanded=True):
        st.selectbox("Observatory", options, key=preset_key,
                     on_change=_apply_preset, args=(key_prefix, preset_key))
        is_manual = st.session_state[preset_key] == MANUAL_OPTION

        c1, c2, c3, c4, c5 = st.columns([2, 1, 1, 1, 1.4])
        name = c1.text_input("Site name", key=f"{key_prefix}_name", disabled=not is_manual)
        lat = c2.number_input("Latitude (deg N)", min_value=-90.0, max_value=90.0,
                               format="%.4f", key=f"{key_prefix}_lat", disabled=not is_manual)
        lon = c3.number_input("Longitude (deg E)", min_value=-180.0, max_value=180.0,
                               format="%.4f", key=f"{key_prefix}_lon", disabled=not is_manual)
        elev = c4.number_input("Elevation (m)", min_value=-500.0,
                                max_value=9000.0, key=f"{key_prefix}_elev", disabled=not is_manual)
        tz = c5.text_input("Timezone (IANA)", key=f"{key_prefix}_tz",
                            disabled=not is_manual,
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
    st.caption("SNR as a function of exposure time, star magnitude, and sky brightness, "
               "for the selected instrument.")

    instrument = instrument_form("snr")

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
                                     exp_fixed=exp_fixed, star_mags3=star_mags, instrument=instrument)
    st.pyplot(fig)


# =============================================================================
# Page 2: Light Curve Simulator
# =============================================================================
LIGHTCURVE_PRESETS = {
    "Custom": None,
    "HAT-P-32 b": dict(star_mag=11.44, k=0.1508, b_imp=0.147, t14_h=3.07, u1=0.30, u2=0.28),
    "HD 189733 b": dict(star_mag=7.65, k=0.1504, b_imp=0.671, t14_h=1.83, u1=0.30, u2=0.28),
    "HD 209458 b": dict(star_mag=7.65, k=0.1209, b_imp=0.507, t14_h=3.07, u1=0.30, u2=0.28),
    # First transiting exoplanet ever discovered (1999); bright host, Dec +18.9.
    "WASP-33 b": dict(star_mag=8.14, k=0.1118, b_imp=0.210, t14_h=2.85, u1=0.30, u2=0.28),
    # Bright, unusual pulsating (delta Scuti) host star; Dec +37.6.
    "WASP-12 b": dict(star_mag=11.57, k=0.1170, b_imp=0.424, t14_h=3.00, u1=0.30, u2=0.28),
    # Famous inspiraling/tidally-disrupting hot Jupiter; Dec +29.7.
    "Qatar-1 b": dict(star_mag=12.69, k=0.1463, b_imp=0.645, t14_h=1.66, u1=0.30, u2=0.28),
    # Deep (~2.1%) transit, circumpolar from mid-northern latitudes at Dec +65.
    "GJ 486 b": dict(star_mag=11.39, k=0.0372, b_imp=0.120, t14_h=1.016, u1=0.30, u2=0.28),
    # Earth-sized (1.3 Re) rocky planet, bright M dwarf host; Dec +9.8 (~50 deg max alt from Frankfurt).
    "K2-18 b": dict(star_mag=13.48, k=0.0538, b_imp=0.65, t14_h=2.97, u1=0.30, u2=0.28),
    # Habitable-zone sub-Neptune (JWST H2O/DMS target); Dec +7.6, well placed from mid-north.
    "LTT 1445 A b": dict(star_mag=11.22, k=0.0451, b_imp=0.17, t14_h=1.37, u1=0.30, u2=0.28),
    # Closest Earth-sized (1.3 Re) transiting rocky planet; Dec -16.3, low but reachable.
    "GJ 1214 b": dict(star_mag=14.71, k=0.1159, b_imp=0.264, t14_h=0.870, u1=0.30, u2=0.28),
    # Classic amateur-reachable sub-Neptune, deep 1.3% transit; Dec +4.9.
    "TRAPPIST-1 b": dict(star_mag=18.80, k=0.0859, b_imp=0.095, t14_h=0.601, u1=0.30, u2=0.28),
    # Innermost of 7 Earth-sized planets; very faint in V (use I band in practice), Dec -5.0.
    "TRAPPIST-1 c": dict(star_mag=18.80, k=0.0844, b_imp=0.109, t14_h=0.701, u1=0.30, u2=0.28),
    # Earth-sized, same ultracool dwarf host as TRAPPIST-1 b; short (42 min) transit.
}


def page_lightcurve_simulator():
    st.header("Light Curve Simulator")
    st.caption("Simulate a single-target transit light curve for given star/planet "
               "parameters and instrument exposure settings.")

    instrument = instrument_form("lc")

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
                                   obs_hours=obs_hours, bin_min=bin_min, seed=int(seed),
                                   instrument=instrument)
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
    instrument = instrument_form("transits")

    c1, c2, c3 = st.columns(3)
    min_alt = c1.slider("Min. altitude throughout transit (deg)", 0.0, 80.0, 35.0, 1.0)
    bright_limit = c1.number_input("Bright limit (V mag, may saturate below)", value=5.5, format="%.1f")
    faint_limit = c2.number_input("Faint limit (V mag, poor SNR above)", value=10.5, format="%.1f")
    min_depth = c2.number_input("Min. transit depth (%)", value=1.5, min_value=0.0, format="%.2f")
    st.session_state.setdefault("transits_sky_mag", site.sky_mag)
    sky_mag_obs = c3.number_input("Sky brightness at your site (mag/arcsec2)",
                                   format="%.2f", key="transits_sky_mag",
                                   help="Prefilled from the preset observatory; edit freely for tonight's actual "
                                        "conditions (moon phase, haze, etc).")
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

        show_weather(dark_window)

        planets = cached_fetch_transiting_planets(bright_limit, faint_limit, min_depth)
        if planets is None:
            st.error("NASA Exoplanet Archive is currently unreachable (the service may be down or "
                     "blocking requests) - please try again later.")
            return
        st.write(f"{len(planets)} confirmed transiting planets in magnitude/depth range retrieved from archive.")

        if include_toi:
            toi_planets = cached_fetch_toi_candidates(bright_limit, faint_limit, min_depth)
            if toi_planets is None:
                st.warning("TESS Objects of Interest (ExoFOP/TOI) table unreachable - "
                           "continuing with confirmed planets only.")
            else:
                st.write(f"{len(toi_planets)} unverified TOI candidates in range retrieved from ExoFOP/TOI table.")
                planets = planets + toi_planets

        top = find_tonight_transits(planets, dark_window, min_alt, sky_mag_obs, n_top, instrument=instrument)

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

        fig = build_transit_figures(top, dark_window, sky_mag_obs, instrument=instrument)
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
    instrument = instrument_form("eb")

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
    st.session_state.setdefault("eb_sky_mag", site.sky_mag)
    sky_mag_obs_eb = c3.number_input("Sky brightness at your site (mag/arcsec2)",
                                      format="%.2f", key="eb_sky_mag",
                                      help="Prefilled from the preset observatory; edit freely for tonight's "
                                           "actual conditions (moon phase, haze, etc).")
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

        show_weather(dark_window)

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
            include_compact_systems=include_compact, instrument=instrument,
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

        fig = build_eb_figures(top, dark_window, sky_mag_obs_eb, instrument=instrument)
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
