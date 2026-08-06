"""
Instrument model (optical tube + camera) and the pure, side-effect-free
photometry/SNR functions built on top of it.

The flux zero point is empirically calibrated for the Seestar S50 specifically
(community reports: a mag-5.25 star saturates its 12-bit sensor in ~1s) and
scaled to other instruments by relative light-collecting power
(aperture^2 x QE) - see `flux_zp_e_s()` - rather than derived from a
from-scratch physical V-band zero point. This keeps Seestar numbers exactly
as calibrated while extending sensibly to larger apertures and different QE.
"""
from dataclasses import dataclass

import numpy as np


@dataclass
class OTA:
    """Optical tube assembly: aperture + focal length set light-collecting
    area and plate scale, independent of whatever camera it's paired with."""
    name: str
    aperture_mm: float
    focal_length_mm: float

    @property
    def f_ratio(self) -> float:
        return self.focal_length_mm / self.aperture_mm


@dataclass
class Camera:
    """Detector: pixel geometry + noise/QE/full-well characteristics."""
    name: str
    pixel_size_um: float
    full_well_e: float      # electrons
    read_noise_e: float
    dark_current_e_s: float
    qe: float                # quantum efficiency, 0-1
    gain_e_per_adu: float
    bit_depth: int = 16


@dataclass
class Instrument:
    ota: OTA
    camera: Camera
    photometric_aperture_radius_arcsec: float = 7.178015   # fixed on-sky aperture radius (= 3 px on the Seestar)

    @property
    def name(self) -> str:
        return f"{self.ota.name} + {self.camera.name}"

    @property
    def pixel_scale_arcsec(self) -> float:
        return (self.camera.pixel_size_um * 1e-3 / self.ota.focal_length_mm) * (180 / np.pi * 3600)

    @property
    def aperture_radius_px(self) -> float:
        return self.photometric_aperture_radius_arcsec / self.pixel_scale_arcsec

    @property
    def n_pix(self) -> float:
        return np.pi * self.aperture_radius_px ** 2


# =============================================================================
# Seestar S50 / Sony IMX462 - the reference instrument
# =============================================================================
SEESTAR_OTA = OTA(name="Seestar S50", aperture_mm=50.0, focal_length_mm=250.0)
SEESTAR_CAMERA = Camera(
    name="Sony IMX462", pixel_size_um=2.9,
    full_well_e=3855 * 1.5,   # SATURATION_ADU (12-bit) * GAIN_E_PER_ADU
    read_noise_e=3.0, dark_current_e_s=0.1, qe=0.8, gain_e_per_adu=1.5, bit_depth=12,
)
SEESTAR = Instrument(ota=SEESTAR_OTA, camera=SEESTAR_CAMERA)

_SEESTAR_SAT_MAG_1S = 5.25   # empirical: mag at which the Seestar saturates in a 1s exposure
_SEESTAR_FLUX_ZP_E_S = SEESTAR_CAMERA.full_well_e / 10 ** (-0.4 * _SEESTAR_SAT_MAG_1S)


def flux_zp_e_s(instrument: Instrument = SEESTAR) -> float:
    """Zero-point flux (e-/s for a mag-0 star), scaled from the Seestar's
    empirical calibration by relative light-collecting power (aperture^2 x QE)."""
    power_ratio = (instrument.ota.aperture_mm ** 2 * instrument.camera.qe) / \
                  (SEESTAR_OTA.aperture_mm ** 2 * SEESTAR_CAMERA.qe)
    return _SEESTAR_FLUX_ZP_E_S * power_ratio


# =============================================================================
# Saturation / exposure
# =============================================================================
def calculate_max_exposure(magnitude, instrument: Instrument = SEESTAR, safety_factor: float = 0.8):
    """Maximum exposure time (s) before a star of `magnitude` reaches
    `safety_factor` of the camera's full well (default: 80%)."""
    target_e = instrument.camera.full_well_e * safety_factor
    return target_e / star_flux(magnitude, instrument)


def calculate_saturation_magnitude(exposure_time_s, instrument: Instrument = SEESTAR,
                                    safety_factor: float = 1.0):
    """Approximate magnitude at which a star reaches `safety_factor` of full
    well (default 100% = true saturation) in a given exposure time."""
    target_e = instrument.camera.full_well_e * safety_factor
    return -2.5 * np.log10(target_e / (flux_zp_e_s(instrument) * exposure_time_s))


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


def star_flux(mag, instrument: Instrument = SEESTAR):
    """Star flux in e-/s."""
    return flux_zp_e_s(instrument) * 10 ** (-0.4 * mag)


def sky_flux_per_pix(sky_mag_arcsec2, instrument: Instrument = SEESTAR):
    """Sky background in e-/pix/s given surface brightness in mag/arcsec^2."""
    return flux_zp_e_s(instrument) * 10 ** (-0.4 * sky_mag_arcsec2) * instrument.pixel_scale_arcsec ** 2


def snr_full(mag, exposure_s, sky_mag_arcsec2=20.5, instrument: Instrument = SEESTAR):
    """Full aperture-photometry SNR including sky, dark, and read noise."""
    S = star_flux(mag, instrument) * exposure_s
    B = sky_flux_per_pix(sky_mag_arcsec2, instrument) * exposure_s
    D = instrument.camera.dark_current_e_s * exposure_s
    noise_sq = S + instrument.n_pix * (B + D + instrument.camera.read_noise_e ** 2)
    return S / np.sqrt(noise_sq)
