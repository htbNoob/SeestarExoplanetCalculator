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
Adjust the **USER PARAMETERS** block in the transit simulation cell to match your target:

| Parameter | Description |
|-----------|-------------|
| `STAR_MAG` | Host star V-band magnitude |
| `K` | Planet-to-star radius ratio (Rp/Rs) |
| `B_IMP` | Impact parameter (0 = central transit) |
| `T14_H` | Total transit duration T₁₄ (hours) |
| `U1`, `U2` | Quadratic limb-darkening coefficients |
| `EXPOSURE_S` | Single-frame exposure time (seconds) |
| `SKY_MAG` | Sky surface brightness (mag/arcsec²) |

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

## Example Output

The notebook includes a worked example for **HAT-P-32 b** (mag 11.3, 2.22% depth), demonstrating that the Seestar S50 can detect the transit with stacked short exposures.

## License

MIT
