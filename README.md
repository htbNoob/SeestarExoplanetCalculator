# SeestarExoplanetCalculator

A Jupyter notebook for planning and simulating exoplanet transit observations with the **ZWO Seestar S50** smart telescope (50 mm aperture, f/5, Sony IMX462 sensor).

## Features

- **Unified configuration cell** — all imports and every tunable parameter in a single cell (cell 2); run it once at startup to configure the entire notebook
- **Saturation calculator** — find the brightest star magnitude a given exposure can handle before saturating the 12-bit sensor
- **Max exposure estimator** — given a target magnitude, compute the longest safe exposure (with configurable safety margin)
- **Full SNR model** — aperture-photometry SNR including shot noise, sky background, dark current, and read noise
- **SNR overview plots** — three diagnostic charts: SNR vs exposure time, SNR vs star magnitude, and SNR vs sky brightness (Bortle scale)
- **Simulated transit light curve** — physically-motivated, limb-darkened transit simulation with per-frame photon noise; produces a publication-style figure with residuals and a stellar disk diagram
- **Tonight's best observable transits** — queries the NASA Exoplanet Archive TAP service live, filters by altitude, magnitude, and depth, ranks candidates, and plots simulated light curves for the top hits
- **Tonight's best eclipsing binaries & compact-object systems** — queries AAVSO VSX for EA/EB/EW variables and checks a curated list of white-dwarf, neutron-star, and black-hole binaries; ranks by depth × SNR × coverage and plots trapezoidal simulated light curves

## Requirements

```
numpy
matplotlib
astropy
requests
```

Install with:

```bash
pip install numpy matplotlib astropy requests
```

## Usage

Open `SeestarExoplanetCalculator.ipynb` in Jupyter or VS Code and **run cell 2 first** — it contains all imports and every tunable parameter. You only need to re-run cell 2 whenever you change a parameter.

### All Parameters (cell 2)

All configuration lives in cell 2. The key sections are:

**Simulated transit light curve (used by cell 8)**

| Parameter | Description |
|-----------|-------------|
| `STAR_MAG` | Host star V-band magnitude |
| `K` | Planet-to-star radius ratio (Rp/Rs) |
| `B_IMP` | Impact parameter (0 = central transit) |
| `T14_H` | Total transit duration T₁₄ (hours) |
| `U1`, `U2` | Quadratic limb-darkening coefficients |
| `EXPOSURE_S` | Single-frame exposure time (seconds) |
| `SKY_MAG` | Sky surface brightness (mag/arcsec²) |

**Tonight's best observable transits (used by cell 10)**

| Parameter | Description |
|-----------|-------------|
| `LAT` / `LON` / `ELEV_M` | Observer coordinates (decimal degrees, metres) |
| `SITE_NAME` | Label for plot titles |
| `MIN_ALT` | Minimum star altitude throughout the transit (°) |
| `BRIGHT_LIMIT` / `FAINT_LIMIT` | V-magnitude window (saturation / SNR limits) |
| `MIN_DEPTH` | Minimum transit depth to consider (%) |
| `SKY_MAG_OBS` | Sky surface brightness — see table below |
| `N_TOP` | How many top-ranked transits to display |

**Tonight's best eclipsing binaries & compact-object systems (used by cell 12)**

| Parameter | Description |
|-----------|-------------|
| `MIN_DEPTH_EB` | Minimum eclipse depth to consider (magnitudes) |
| `MAX_PERIOD_EB` | Maximum orbital period (days) |
| `MIN_ALT_EB` | Minimum star altitude throughout eclipse (°) |
| `N_TOP_EB` | How many top-ranked systems to display |
| `BIN_MIN_EB` | Bin width for the simulated light curve (minutes) |

**Sky brightness guide (`SKY_MAG_OBS` / `SKY_MAG`):**

| Value | Bortle class | Typical location |
|-------|-------------|-----------------|
| 17.0 | 9 | Inner city |
| 18.5 | 7 | Suburban |
| 20.5 | 4–5 | Rural / suburban border |
| 21.7 | 1–2 | Truly dark site |

## Instrument Constants (Seestar S50 / IMX462)

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

## Cells at a Glance

| # | Title | What it does |
|---|-------|-------------|
| 1 | Title | Notebook header (markdown) |
| 2 | **Imports & Parameters** | All imports + every tunable parameter — **run this first** |
| 3 | Saturation & Photometry Calculator | Header (markdown) |
| 4 | Core Photometry Functions | Saturation magnitude, max exposure, SNR, transit depth conversions |
| 5 | SNR Model | Header with formula (markdown) |
| 6 | SNR Overview Plots | Three diagnostic charts: SNR vs exposure time, SNR vs star magnitude, SNR vs sky brightness |
| 7 | Simulated Transit Light Curve | Header (markdown) |
| 8 | Simulated Transit Light Curve | Full limb-darkened transit simulation with photon noise, residuals panel, and stellar disk diagram |
| 9 | Tonight's Best Observable Transits | Header (markdown) |
| 10 | **Tonight's Best Observable Transits** | Queries NASA Exoplanet Archive for tonight's transits at your location, ranks them, and plots simulated light curves |
| 11 | Tonight's Best Eclipsing Binaries | Header (markdown) |
| 12 | **Tonight's Best Eclipsing Binaries & Compact-Object Systems** | Queries AAVSO VSX for EA/EB/EW stars + curated WD/NS/BH binaries; ranks and plots trapezoidal light curves |

## Example Output

The notebook includes a worked example for **HAT-P-32 b** (mag 11.3, 2.22% depth), demonstrating that the Seestar S50 can detect the transit with stacked short exposures.

The **Tonight's Transits** cell (cell 10) fetches live data from the [NASA Exoplanet Archive TAP service](https://exoplanetarchive.ipac.caltech.edu/TAP/sync) (no API key needed) and filters by altitude, magnitude range, and transit depth for your site.

The **Tonight's Eclipsing Binaries** cell (cell 12) queries the [AAVSO VSX API](https://www.aavso.org/vsx/) for eclipsing variables and always includes a curated set of compact-object binaries (Cyg X-1, HZ Her, V471 Tau, HW Vir, and others) regardless of live query results.

## License

MIT
