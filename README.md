# Seestar Exoplanet Transit Observation Planner

A Streamlit web app for planning and simulating exoplanet transit (and eclipsing-binary) observations with the **ZWO Seestar S50** smart telescope (50 mm aperture, f/5, Sony IMX462 sensor).

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
SNR as a function of exposure time, star magnitude, and sky brightness, for the Seestar S50's aperture-photometry model (shot noise, sky background, dark current, read noise). Three linked panels, each with adjustable overlay values.

### 2. Light Curve Simulator
Simulates a single-target transit light curve: limb-darkened stellar disk, per-frame photon noise, binned overlay, residuals panel, and a stellar-disk diagram showing the planet's path. Inputs: star magnitude, Rp/Rs, impact parameter, T₁₄, limb-darkening coefficients, exposure time, sky brightness, observation window, bin width, and random seed. Includes presets for HAT-P-32 b and HD 189733 b.

### 3. Tonight's Best Observable Transits
Queries the [NASA Exoplanet Archive TAP service](https://exoplanetarchive.ipac.caltech.edu/TAP/sync) live (no API key needed), filters by your site's altitude/magnitude/depth constraints, and ranks candidates by depth × SNR × sky coverage × altitude. Plots simulated light curves for the top hits.

Optionally include **unverified TOI candidates** (TESS Objects of Interest — the same catalog ExoFOP-TESS serves) alongside confirmed planets via a checkbox. These are flagged separately in the results table and plot titles, since their transit parameters are less certain and some will turn out to be false positives.

### 4. Tonight's Best Eclipsing Binaries & Compact-Object Systems
Queries [AAVSO VSX](https://www.aavso.org/vsx/) for eclipsing-type variables (EA/EB/EW) in a cone around your zenith at local midnight, plus a curated list of white-dwarf, neutron-star, and black-hole binaries (V471 Tau, HW Vir, HZ Her, Cyg X-1, and others) that's included regardless of live query results. Ranks and plots trapezoidal eclipse models for the top systems.

**Note:** the VSX query can take 1–2 minutes on dense fields — this is a known server-side limitation (see comments in `seestar_eb.py`), not an app bug.

### Observer site & local time
Pages 3 & 4 share a site form (name, latitude, longitude, elevation, IANA timezone e.g. `Europe/Berlin`) persisted across both pages. "Tonight's dark window" is computed as the time the sun is below −12° (nautical twilight) at your site, shown in both UTC and your local timezone.

## Instrument constants (Seestar S50 / IMX462)

| Parameter | Value |
|-----------|-------|
| Aperture | 50 mm |
| Focal length | 250 mm (f/5) |
| Pixel size | 2.9 µm |
| Pixel scale | ~2.39 arcsec/pixel |
| Saturation | 3855 ADU (12-bit) |
| Read noise | ~3 e⁻ RMS |
| Gain | ~1.5 e⁻/ADU |
| QE | ~80% (STARVIS estimate) |

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
| `seestar_core.py` | Instrument constants + pure physics functions (saturation, SNR, mag↔depth conversions) — no I/O |
| `seestar_site.py` | Observer site, tonight's dark-window computation, UTC↔local time formatting |
| `seestar_snr.py` | SNR overview plot |
| `seestar_lightcurve.py` | Single-target transit simulator |
| `seestar_transits.py` | NASA Exoplanet Archive query (confirmed planets + TOI candidates), ranking, and plotting |
| `seestar_eb.py` | AAVSO VSX query, curated compact-object list, ranking, and plotting |

## License

MIT
