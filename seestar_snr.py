"""SNR overview plot (SNR vs exposure / magnitude / sky background)."""
import matplotlib.pyplot as plt
import numpy as np

from seestar_core import calculate_max_exposure, calculate_saturation_magnitude, snr_full

BORTLE_LABELS = {17.0: "Bortle 9\n(inner city)", 18.5: "Bortle 7",
                  20.5: "Bortle 4-5", 21.7: "Bortle 1-2\n(dark)"}


def build_snr_overview_figure(ref_sky=20.5, star_mags=(8, 10, 12, 13, 14, 15),
                               exp_times=(1, 5, 10, 30, 60, 120),
                               exp_fixed=10, star_mags3=(8, 10, 12, 13, 14, 15)):
    """Reproduces the notebook's 3-panel SNR overview, parameterized for the UI."""
    fig, axes = plt.subplots(1, 3, figsize=(16, 5))
    fig.suptitle("Seestar S50 - SNR Overview (IMX462, 50 mm, f/5)", fontsize=13, fontweight="bold")

    colors = plt.cm.plasma_r(np.linspace(0.15, 0.85, max(len(star_mags), 1)))

    # Panel 1: SNR vs exposure time
    ax1 = axes[0]
    t_arr = np.logspace(-1, 3, 400)
    for i, mag in enumerate(star_mags):
        snr_arr = snr_full(mag, t_arr, sky_mag_arcsec2=ref_sky)
        sat_t = calculate_max_exposure(mag)
        mask = t_arr <= sat_t
        ax1.loglog(t_arr[mask], snr_arr[mask], color=colors[i], lw=2, label=f"mag {mag}")
        if np.any(~mask):
            ax1.axvline(sat_t, color=colors[i], lw=0.7, ls=":", alpha=0.6)
    ax1.axhline(100, color="gray", ls="--", lw=1, alpha=0.7, label="SNR = 100")
    ax1.set_xlabel("Exposure time (s)")
    ax1.set_ylabel("SNR")
    ax1.set_title(f"SNR vs Exposure Time\n(sky = {ref_sky:.1f} mag/arcsec2)")
    ax1.legend(fontsize=8, loc="upper left")
    ax1.grid(True, which="both", alpha=0.2)
    ax1.set_xlim(0.1, 1000)

    # Panel 2: SNR vs star magnitude
    ax2 = axes[1]
    mag_arr = np.linspace(5, 16, 300)
    colors2 = plt.cm.viridis(np.linspace(0.1, 0.9, max(len(exp_times), 1)))
    for i, exp in enumerate(exp_times):
        snr_arr = snr_full(mag_arr, exp, sky_mag_arcsec2=ref_sky)
        sat_mag_limit = calculate_saturation_magnitude(exp)
        mask = mag_arr >= sat_mag_limit
        ax2.semilogy(mag_arr[mask], snr_arr[mask], color=colors2[i], lw=2, label=f"{exp} s")
    ax2.axhline(100, color="gray", ls="--", lw=1, alpha=0.7, label="SNR = 100")
    ax2.set_xlabel("Star magnitude (Vega)")
    ax2.set_ylabel("SNR")
    ax2.set_title(f"SNR vs Star Magnitude\n(sky = {ref_sky:.1f} mag/arcsec2)")
    ax2.legend(fontsize=8, title="Exposure", loc="upper right")
    ax2.grid(True, which="both", alpha=0.2)
    ax2.set_xlim(5, 16)

    # Panel 3: SNR vs sky background
    ax3 = axes[2]
    sky_arr = np.linspace(16, 22, 300)
    colors3 = plt.cm.plasma_r(np.linspace(0.15, 0.85, max(len(star_mags3), 1)))
    for sky_b, label in BORTLE_LABELS.items():
        ax3.axvline(sky_b, color="lightgray", lw=0.8, ls=":")
        ax3.text(sky_b + 0.05, 1.5, label, fontsize=6, color="gray", va="bottom")
    for i, mag in enumerate(star_mags3):
        snr_arr = np.array([snr_full(mag, exp_fixed, sky_mag_arcsec2=s) for s in sky_arr])
        ax3.semilogy(sky_arr, snr_arr, color=colors3[i], lw=2, label=f"mag {mag}")
    ax3.axhline(100, color="gray", ls="--", lw=1, alpha=0.7, label="SNR = 100")
    ax3.set_xlabel("Sky brightness (mag/arcsec2)  ->  darker")
    ax3.set_ylabel("SNR")
    ax3.set_title(f"SNR vs Sky Background\n(t = {exp_fixed:.0f} s)")
    ax3.legend(fontsize=8, loc="upper left")
    ax3.grid(True, which="both", alpha=0.2)
    ax3.set_xlim(16, 22)

    plt.tight_layout()
    return fig
