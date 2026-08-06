# Seestar Exoplanet Transit Observation Planner

A Streamlit web app for planning and simulating exoplanet transit (and eclipsing-binary) observations — built around the **ZWO Seestar S50** but extensible to other telescope + camera combinations (see "Instrument model" below).

Four interactive tools in one app (`app.py`), backed by a set of plain-Python calculation modules (`seestar_*.py`).

## Requirements

```
streamlit
numpy
matplotlib
astropy
requests
```

Install with:

```bash
pip install -r requirements.txt
```

## Running the app

```bash
streamlit run app.py
```

This opens the app in your browser (default `http://localhost:8501`). Use the sidebar to switch between the four tools below. Stop it with `Ctrl+C`.

## The four tools

### 1. SNR Explorer
SNR as a function of exposure time, star magnitude, and sky brightness, for the selected instrument's aperture-photometry model (shot noise, sky background, dark current, read noise). Three linked panels, each with adjustable overlay values.

### 2. Light Curve Simulator
Simulates a single-target transit light curve: limb-darkened stellar disk, per-frame photon noise, binned overlay, residuals panel, and a stellar-disk diagram showing the planet's path. Inputs: star magnitude, Rp/Rs, impact parameter, T₁₄, limb-darkening coefficients, exposure time, sky brightness, observation window, bin width, and random seed.

Includes presets for six well-established, northern-hemisphere-visible targets (parameters pulled live from the NASA Exoplanet Archive): HAT-P-32 b, HD 189733 b, HD 209458 b (the first transiting exoplanet ever discovered), WASP-33 b (bright, pulsating host star), WASP-12 b (famous inspiraling hot Jupiter), and Qatar-1 b (circumpolar from mid-northern latitudes, Dec +65°).

### 3. Tonight's Best Observable Transits
Queries the [NASA Exoplanet Archive TAP service](https://exoplanetarchive.ipac.caltech.edu/TAP/sync) live (no API key needed), filters by your site's altitude/magnitude/depth constraints, and ranks candidates by depth × SNR × sky coverage × altitude. Plots simulated light curves for the top hits.

Optionally include **unverified TOI candidates** (TESS Objects of Interest — the same catalog ExoFOP-TESS serves) alongside confirmed planets via a checkbox. These are flagged separately in the results table and plot titles, since their transit parameters are less certain and some will turn out to be false positives.

### 4. Tonight's Best Eclipsing Binaries & Compact-Object Systems
Queries [AAVSO VSX](https://www.aavso.org/vsx/) for eclipsing-type variables (EA/EB/EW) in a cone around your zenith at local midnight, plus a curated list of white-dwarf, neutron-star, and black-hole binaries (V471 Tau, HW Vir, HZ Her, Cyg X-1, and others) that's included regardless of live query results. Ranks and plots trapezoidal eclipse models for the top systems.

**Note:** the VSX query can take 1–2 minutes on dense fields — this is a known server-side limitation (see comments in `seestar_eb.py`), not an app bug.

### Observer site & local time
Pages 3 & 4 share a site form (name, latitude, longitude, elevation, IANA timezone e.g. `Europe/Berlin`) persisted across both pages. Pick a preset observatory from the dropdown, or choose **Manual (custom location)** to type in your own coordinates:

| Preset | Lat, Lon | Elevation | Timezone | Sky brightness |
|--------|----------|-----------|----------|----------------|
| Frankfurt am Main, Germany | 50.0853, 8.5822 | 100 m | Europe/Berlin | 18.0 mag/arcsec² |
| Taunus Observatory (Kleiner Feldberg) | 50.2217, 8.4458 | 826 m | Europe/Berlin | 20.8 mag/arcsec² (Bortle 4) |
| Teide Observatory, Tenerife | 28.3000, −16.5097 | 2390 m | Atlantic/Canary | 21.59 mag/arcsec² (Bortle 4) |
| Madrona Peak Observatory | 29.8442, −99.3558 | 379 m | America/Chicago | 21.89 mag/arcsec² (Bortle 2) |

Selecting a preset also prefills the "sky brightness" field on the transits/EB pages with that site's typical value — it stays editable, since actual conditions (moon phase, haze) vary night to night.

"Tonight's dark window" is computed as the time the sun is below −12° (nautical twilight) at your site, shown in both UTC and your local timezone.

### Weather forecast
Both "tonight" pages fetch an hourly cloud cover / precipitation probability / temperature / dew point / wind speed forecast for your site from [Open-Meteo](https://open-meteo.com) (free, no API key) and summarize it over just the dark-window hours: a Good/Fair/Poor rating (based on mean cloud cover and precipitation probability) plus a two-panel forecast chart with the dark window shaded. Cached for 30 minutes per site.

Light pollution is still a manual input (the "sky brightness" field) — there's no equivalent free live API for that; it would require bundling a static VIIRS-derived raster dataset, which is a bigger addition than this app currently needs.

## Instrument model

Every page has an **Instrument** picker: pick an optical tube (aperture + focal length) and a camera (pixel size, full well, read noise, dark current, QE) independently, mix and match. The Seestar S50 is the reference/default and its numbers are unchanged from the original calculator; other instruments' flux zero point is derived by scaling the Seestar's empirical calibration by relative light-collecting power (aperture² × QE) — see `flux_zp_e_s()` in `seestar_core.py`.

**Optical tubes:**

| Preset | Aperture | Focal length | f-ratio |
|--------|----------|--------------|---------|
| Seestar S50 | 50 mm | 250 mm | f/5 |
| Taunus Westkuppel (PlaneWave CDK12.5) | 317.5 mm | 2531 mm | f/8 |
| Taunus Ostkuppel 60cm, primary focus | 600 mm | 2000 mm | f/3.3 |
| Taunus Ostkuppel 60cm, secondary focus | 600 mm | 6000 mm | f/10 |
| MPO61 RCOS 610 (24") | 609.6 mm | 4814 mm | f/7.9 |

**Cameras:**

| Preset | Pixel size | Full well | Read noise | QE | Notes |
|--------|-----------|-----------|------------|-----|-------|
| Seestar Sony IMX462 | 2.9 µm | 5783 e⁻ | 3.0 e⁻ | 80% | Reference calibration |
| Atik 383L+ (KAF-8300 CCD) | 5.4 µm | 25500 e⁻ | 8.5 e⁻ | 50% | Dark current & gain estimated |
| Moravian C3-61000 (Sony IMX455 CMOS) | 3.76 µm | 50000 e⁻ | 2.0 e⁻ | 85% | Dark current & gain estimated |
| QHY461PH (MPO61) | 3.76 µm | 75600 e⁻ | 3.5 e⁻ (est.) | 85% (est.) | Full well & gain from MPO61's own config; read noise/QE/dark current estimated from QHYCCD literature |

Values not published directly in the source spec sheets (dark current for both Taunus cameras; gain for Atik/Moravian; read noise/QE/dark current for the QHY461PH) are estimated and flagged in `app.py`'s `CAMERA_PRESETS` comments — worth confirming against your own calibration frames if precision matters.

**Not modeled:** per-filter photometry (all instruments share the single V-band-proxy zero point the Seestar was calibrated on — swapping the camera/aperture changes signal levels correctly, but not filter response), and sensor binning (e.g. MPO61 commonly runs 2×2 binned — pick the native pixel size and mentally adjust, or ask for binning support as a follow-up).

## Sky brightness guide

Used for the "sky brightness" inputs throughout the app:

| Value (mag/arcsec²) | Bortle class | Typical location |
|----------------------|---------------|-------------------|
| 17.0 | 9 | Inner city |
| 18.5 | 7 | Suburban |
| 20.5 | 4–5 | Rural / suburban border |
| 21.7 | 1–2 | Truly dark site |

## Code layout

| File | Contents |
|------|----------|
| `app.py` | Streamlit UI — page routing, forms, cached network calls |
| `seestar_core.py` | OTA/Camera/Instrument model + pure physics functions (saturation, SNR, mag↔depth conversions) — no I/O |
| `seestar_site.py` | Observer site, tonight's dark-window computation, UTC↔local time formatting |
| `seestar_snr.py` | SNR overview plot |
| `seestar_lightcurve.py` | Single-target transit simulator |
| `seestar_transits.py` | NASA Exoplanet Archive query (confirmed planets + TOI candidates), ranking, and plotting |
| `seestar_eb.py` | AAVSO VSX query, curated compact-object list, ranking, and plotting |
| `seestar_weather.py` | Open-Meteo forecast fetch, dark-window summary/rating, forecast plot |

## License

MIT
