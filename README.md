# SeestarExoplanetCalculator

A Jupyter notebook for planning and simulating exoplanet transit observations with the **ZWO Seestar S50** smart telescope (50 mm aperture, f/5, Sony IMX462 sensor).

## Features

- **Saturation calculator** — find the brightest star magnitude a given exposure can handle before saturating the 12-bit sensor
- **Max exposure estimator** — given a target magnitude, compute the longest safe exposure (with configurable safety margin)
- **Full SNR model** — aperture-photometry SNR including shot noise, sky background, dark current, and read noise
- **SNR overview plots** — three diagnostic charts: SNR vs exposure time, SNR vs star magnitude, and SNR vs sky brightness (Bortle scale)
- **Simulated transit light curve** — physically-motivated, limb-darkened transit simulation with per-frame photon noise; produces a publication-style figure with residuals and a stellar disk diagram

## Requirements

```
numpy
matplotlib
```

Install with:

```bash
pip install numpy matplotlib
```

## Usage

Open `SeestarExoplanetCalculator.ipynb` in Jupyter or VS Code and run the cells in order.

### Simulated Transit Light Curve (cell 5)

Adjust the **USER PARAMETERS** block to match your target:

| Parameter | Description |
|-----------|-------------|
| `STAR_MAG` | Host star V-band magnitude |
| `K` | Planet-to-star radius ratio (Rp/Rs) |
| `B_IMP` | Impact parameter (0 = central transit) |
| `T14_H` | Total transit duration T₁₄ (hours) |
| `U1`, `U2` | Quadratic limb-darkening coefficients |
| `EXPOSURE_S` | Single-frame exposure time (seconds) |
| `SKY_MAG` | Sky surface brightness (mag/arcsec²) |

### Tonight's Best Observable Transits (cell 7)

Set your location and site conditions at the top of the cell:

| Parameter | Description |
|-----------|-------------|
| `LAT` / `LON` / `ELEV_M` | Observer coordinates (decimal degrees, metres) |
| `SITE_NAME` | Label for plot titles |
| `MIN_ALT` | Minimum star altitude throughout the transit (°) |
| `BRIGHT_LIMIT` / `FAINT_LIMIT` | V-magnitude window (saturation / SNR limits) |
| `MIN_DEPTH` | Minimum transit depth to consider (%) |
| `SKY_MAG_OBS` | Sky surface brightness — see table below |
| `N_TOP` | How many top-ranked transits to display |

**Sky brightness guide (`SKY_MAG_OBS`):**

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
| 2 | Saturation & Photometry Calculator | Core functions: saturation magnitude, max exposure, SNR, transit depth conversions |
| 4 | SNR Overview Plots | Three diagnostic charts: SNR vs exposure time, SNR vs star magnitude, SNR vs sky brightness |
| 5 | Simulated Transit Light Curve | Full limb-darkened transit simulation with photon noise, residuals panel, and stellar disk diagram |
| 7 | **Tonight's Best Observable Transits** | Queries NASA Exoplanet Archive for tonight's transits at your location, ranks them, and plots simulated light curves |

## Example Output

The notebook includes a worked example for **HAT-P-32 b** (mag 11.3, 2.22% depth), demonstrating that the Seestar S50 can detect the transit with stacked short exposures.

The **Tonight's Transits** cell fetches live data from the [NASA Exoplanet Archive TAP service](https://exoplanetarchive.ipac.caltech.edu/TAP/sync) (no API key needed) and filters by altitude, magnitude range, and transit depth for your site.

## License

MIT
