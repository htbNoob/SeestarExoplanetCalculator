"""
Core Seestar S50 / Sony IMX462 instrument model.

Pure, side-effect-free functions and constants shared by the notebook and the
Streamlit app. Ported from the calculator notebook's setup/SNR cells without
behavioural changes.
"""
import numpy as np

# =============================================================================
# SEESTAR S50 / SONY IMX462 instrument constants
# =============================================================================
APERTURE_MM = 50.0          # mm
FOCAL_LENGTH_MM = 250.0     # mm
F_RATIO = 5.0                # f/5
PIXEL_SIZE_UM = 2.9          # um
SENSOR_AREA_MM2 = np.pi * (APERTURE_MM / 2) ** 2  # mm^2
PIXEL_AREA_M2 = (PIXEL_SIZE_UM * 1e-6) ** 2         # m^2
SATURATION_ADU = 3855        # 12-bit saturation (Sony datasheet)
QE = 0.8                     # Quantum efficiency (STARVIS estimate)
VEGA_ZERO_POINT = 3631       # Photons/cm^2/s/A for V-band (Vega magnitude 0)

GAIN_E_PER_ADU = 1.5         # e-/ADU (typical, from Seestar FITS headers)

# Pixel scale in arcsec/pixel
PIXEL_SCALE_ARCSEC = (PIXEL_SIZE_UM * 1e-3 / FOCAL_LENGTH_MM) * (180 / np.pi * 3600)

READ_NOISE_E = 3.0        # e- RMS (IMX462 spec)
DARK_CURRENT_E_S = 0.1    # e-/pix/s

# Photometric aperture: radius = 3 pix (~7 arcsec), typical for this pixel scale
APERTURE_RADIUS_PIX = 3.0
N_PIX = np.pi * APERTURE_RADIUS_PIX ** 2   # ~28 pixels

# Empirical zero-point: at 1s, a mag-5.25 star fills the well
FWC_E = SATURATION_ADU * GAIN_E_PER_ADU                    # ~5783 e-
FLUX_ZP_E_S = FWC_E / 10 ** (-0.4 * 5.25)                    # e-/s for mag 0 (~728 ke-/s)


# =============================================================================
# Saturation / exposure
# =============================================================================
def calculate_saturation_magnitude(exposure_time_s, gain=GAIN_E_PER_ADU,
                                    aperture_mm=APERTURE_MM,
                                    light_pollution_offset=0.0):
    """Approximate magnitude at which a star saturates the sensor."""
    base_sat_mag_1s = 5.25
    sat_mag = base_sat_mag_1s + 2.5 * np.log10(exposure_time_s)
    sat_mag += light_pollution_offset
    return sat_mag


def calculate_max_exposure(magnitude, gain=GAIN_E_PER_ADU,
                            aperture_mm=APERTURE_MM,
                            light_pollution_offset=0.0,
                            safety_factor=0.8):
    """Maximum exposure time (s) before saturation for a given magnitude."""
    base_sat_mag_1s = 5.25 + light_pollution_offset
    magnitude_adjusted = magnitude + 2.5 * np.log10(1.0 / safety_factor)
    max_exposure = 10 ** ((magnitude_adjusted - base_sat_mag_1s) / 2.5)
    return max_exposure


# =============================================================================
# Magnitude <-> depth conversions
# =============================================================================
def transit_depth_to_magnitude(depth_frac):
    """Fractional transit depth -> magnitude change (negative = dimming)."""
    return -2.5 * np.log10(1 - depth_frac)


def magnitude_to_transit_depth(mag_change):
    """Magnitude change -> fractional depth (0 to 1)."""
    return 1 - 10 ** (-mag_change / 2.5)


# =============================================================================
# SNR
# =============================================================================
def signal_to_noise(flux_electrons, read_noise=3.0, dark_current=0.1):
    """SNR for a given signal in electrons (simple point-source model)."""
    noise_electrons = np.sqrt(flux_electrons + read_noise ** 2)
    return flux_electrons / noise_electrons


def required_flux_for_snr(snr_target, read_noise=3.0):
    """Signal (electrons) required to reach a target SNR."""
    a = 1
    b = -snr_target ** 2
    c = -snr_target ** 2 * read_noise ** 2
    flux = (-b + np.sqrt(b ** 2 - 4 * a * c)) / (2 * a)
    return flux


def star_flux(mag):
    """Star flux in e-/s."""
    return FLUX_ZP_E_S * 10 ** (-0.4 * mag)


def sky_flux_per_pix(sky_mag_arcsec2):
    """Sky background in e-/pix/s given surface brightness in mag/arcsec^2."""
    return FLUX_ZP_E_S * 10 ** (-0.4 * sky_mag_arcsec2) * PIXEL_SCALE_ARCSEC ** 2


def snr_full(mag, exposure_s, sky_mag_arcsec2=20.5,
             read_noise=READ_NOISE_E, dark=DARK_CURRENT_E_S, n_pix=N_PIX):
    """Full aperture-photometry SNR including sky, dark, and read noise."""
    S = star_flux(mag) * exposure_s
    B = sky_flux_per_pix(sky_mag_arcsec2) * exposure_s
    D = dark * exposure_s
    noise_sq = S + n_pix * (B + D + read_noise ** 2)
    return S / np.sqrt(noise_sq)
