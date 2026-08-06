"""Single-target simulated transit light curve (limb-darkened disk model)."""
from dataclasses import dataclass

import matplotlib.gridspec as gridspec
import matplotlib.pyplot as plt
import numpy as np

from seestar_core import calculate_max_exposure, snr_full, transit_depth_to_magnitude


@dataclass
class LightCurveResult:
    fig: plt.Figure
    depth_pct: float
    depth_ppm: float
    depth_mmag: float
    t14_h: float
    t23_h: float
    ingress_min: float
    snr_per_frame: float
    noise_per_frame_mmag: float
    depth_over_sigma: float
    max_safe_exposure_s: float
    frames_in_transit: int
    bins_in_transit: int
    binned_sigma_mmag: float
    grazing: bool


def simulate_transit(star_mag, k, b_imp, t14_h, u1, u2, exposure_s, sky_mag,
                      obs_hours, bin_min, seed=42) -> LightCurveResult:
    """Simulate and plot a transit light curve for a single target.

    Parameters mirror the notebook's cell-7 globals (STAR_MAG, K, B_IMP,
    T14_H, U1, U2, EXPOSURE_S, SKY_MAG, OBS_HOURS, BIN_MIN, SEED).
    """
    # ---- Derived geometry ----
    depth = k ** 2
    denom = np.sqrt((1 + k) ** 2 - b_imp ** 2)
    t23_h = t14_h * np.sqrt(max(0.0, (1 - k) ** 2 - b_imp ** 2)) / denom
    ingress_m = (t14_h - t23_h) / 2 * 60

    snr_frame = snr_full(star_mag, exposure_s, sky_mag_arcsec2=sky_mag)
    sigma_flux = 1.0 / snr_frame
    depth_mmag = transit_depth_to_magnitude(depth) * 1000

    n_in = int(t14_h * 3600 / exposure_s)
    n_bin_in = max(1, int(t14_h * 60 / bin_min))
    binned_sigma_mmag = sigma_flux / np.sqrt(max(n_in, 1) / n_bin_in) * 1e3

    # ---- Limb-darkened stellar disk (computed once) ----
    ng = 200
    xg = np.linspace(-1.0, 1.0, ng)
    Xg, Yg = np.meshgrid(xg, xg)
    Rg = np.sqrt(Xg ** 2 + Yg ** 2)
    mu_g = np.where(Rg <= 1.0, np.sqrt(np.clip(1.0 - Rg ** 2, 0, None)), 0.0)
    I_disk = np.where(Rg <= 1.0, 1.0 - u1 * (1 - mu_g) - u2 * (1 - mu_g) ** 2, 0.0)
    I_total = I_disk.sum()

    # ---- Time array and separation z(t) ----
    n_frames = int(obs_hours * 3600 / exposure_s) + 1
    t_h = np.linspace(-obs_hours / 2, obs_hours / 2, n_frames)
    t_m = t_h * 60.0
    z_arr = np.sqrt(b_imp ** 2 + (2.0 * t_h / t14_h * denom) ** 2)

    # ---- Transit flux (batched) ----
    batch = 50
    flux_model = np.ones(n_frames)
    for i0 in range(0, n_frames, batch):
        zb = z_arr[i0:i0 + batch]
        if np.all(zb >= 1.0 + k):
            continue
        d2 = (Xg[:, :, None] - zb[None, None, :]) ** 2 + Yg[:, :, None] ** 2
        blk = np.sum(I_disk[:, :, None] * (d2 <= k ** 2), axis=(0, 1))
        flux_model[i0:i0 + batch] = 1.0 - blk / I_total

    # ---- Photon noise ----
    rng = np.random.default_rng(seed)
    flux_obs = flux_model + rng.normal(0.0, sigma_flux, n_frames)

    # ---- Bin the light curve ----
    bin_edges = np.arange(t_m[0], t_m[-1] + bin_min, bin_min)
    bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
    bin_flux, bin_err = [], []
    for lo, hi in zip(bin_edges[:-1], bin_edges[1:]):
        idx = (t_m >= lo) & (t_m < hi)
        n = idx.sum()
        if n > 0:
            bin_flux.append(flux_obs[idx].mean())
            bin_err.append(sigma_flux / np.sqrt(n))
        else:
            bin_flux.append(np.nan)
            bin_err.append(np.nan)
    bin_flux = np.array(bin_flux)
    bin_err = np.array(bin_err)
    bin_model = np.interp(bin_centers / 60, t_h, flux_model)

    # ---- Figure ----
    fig = plt.figure(figsize=(16, 8))
    fig.suptitle(
        f"Simulated Transit - Seestar S50  |  "
        f"Star mag {star_mag},  K = {k:.3f} ({depth * 100:.2f}% depth),  "
        f"b = {b_imp},  T14 = {t14_h:.1f} h,  "
        f"t_exp = {exposure_s:.0f} s,  sky = {sky_mag} mag/arcsec2",
        fontsize=10)

    gs = gridspec.GridSpec(2, 3, height_ratios=[3.5, 1.5], hspace=0.06, wspace=0.32)
    ax_lc = fig.add_subplot(gs[0, 0:2])
    ax_res = fig.add_subplot(gs[1, 0:2], sharex=ax_lc)
    ax_disk = fig.add_subplot(gs[:, 2], aspect="equal")
    plt.setp(ax_lc.get_xticklabels(), visible=False)

    ax_lc.errorbar(t_m, flux_obs, yerr=sigma_flux, fmt=".", ms=1.5,
                    color="0.65", elinewidth=0.3, alpha=0.30, label="Single frames", zorder=1)
    ax_lc.errorbar(bin_centers, bin_flux, yerr=bin_err, fmt="o", ms=5,
                    color="steelblue", capsize=3, elinewidth=1.2,
                    label=f"{bin_min:.0f}-min bins", zorder=3)
    ax_lc.plot(t_m, flux_model, "r-", lw=2.0, zorder=4, label="Model (noiseless)")

    y_label = 1 - depth * 1.7
    for tc, lab in [(-t14_h / 2 * 60, "T1"), (-t23_h / 2 * 60, "T2"),
                     (t23_h / 2 * 60, "T3"), (t14_h / 2 * 60, "T4")]:
        ax_lc.axvline(tc, color="tomato", lw=0.9, ls="--", alpha=0.65)
        ax_lc.text(tc, y_label, lab, ha="center", va="top", fontsize=8.5, color="tomato")

    ax_lc.set_ylabel("Relative flux", fontsize=10)
    ax_lc.set_title("Transit Light Curve", fontsize=11)
    ax_lc.legend(fontsize=9, loc="lower center", ncol=3)
    ax_lc.grid(alpha=0.2)
    ax_lc.set_xlim(t_m[0], t_m[-1])

    resid_raw = (flux_obs - flux_model) * 1e3
    resid_bin = (bin_flux - bin_model) * 1e3
    ax_res.errorbar(t_m, resid_raw, yerr=sigma_flux * 1e3, fmt=".", ms=1.5,
                     color="0.65", elinewidth=0.3, alpha=0.30, zorder=1)
    ax_res.errorbar(bin_centers, resid_bin, yerr=bin_err * 1e3, fmt="o", ms=5,
                     color="steelblue", capsize=3, elinewidth=1.2, zorder=3)
    ax_res.axhline(0, color="r", lw=1.5)
    for sign in (+1, -1):
        ax_res.axhline(sign * sigma_flux * 1e3, color="C0", ls="--", lw=0.9, alpha=0.6)

    oot_mask = np.abs(t_h) > t14_h / 2
    rms_oot = np.std(resid_raw[oot_mask]) if oot_mask.any() else np.nan
    ax_res.set_xlabel("Time from mid-transit (min)", fontsize=10)
    ax_res.set_ylabel("Residuals\n(mmag)", fontsize=9)
    ax_res.set_title(
        f"Residuals  -  OOT RMS = {rms_oot:.2f} mmag  "
        f"(expected 1sigma = {sigma_flux * 1e3:.2f} mmag per frame)", fontsize=9)
    ax_res.grid(alpha=0.2)

    ax_disk.imshow(I_disk, extent=[-1, 1, -1, 1], origin="lower",
                    cmap="afmhot", vmin=0, interpolation="bilinear")
    theta = np.linspace(0, 2 * np.pi, 400)
    ax_disk.plot(np.cos(theta), np.sin(theta), "w-", lw=0.9, alpha=0.35)
    x_tr = np.linspace(-denom, denom, 300)
    ax_disk.plot(x_tr, np.full_like(x_tr, b_imp), "w--", lw=1.1, alpha=0.75, label="Planet path")
    ax_disk.add_patch(plt.Circle((0, b_imp), k, color="deepskyblue",
                                  alpha=0.65, lw=1.5, label=f"Planet (K={k})"))
    for x_c in [-denom, denom]:
        ax_disk.add_patch(plt.Circle((x_c, b_imp), k, color="deepskyblue",
                                      alpha=0.20, lw=1.0, ls="--", fill=False))
    ax_disk.set_xlim(-1.3, 1.3)
    ax_disk.set_ylim(-1.3, 1.3)
    ax_disk.set_xlabel("x / R*", fontsize=9)
    ax_disk.set_ylabel("y / R*", fontsize=9)
    ax_disk.set_title(f"Limb-darkened disk\n(u1={u1}, u2={u2})", fontsize=9)
    ax_disk.legend(fontsize=7.5, loc="lower right")
    ax_disk.tick_params(labelsize=8)

    plt.tight_layout()

    return LightCurveResult(
        fig=fig,
        depth_pct=depth * 100,
        depth_ppm=depth * 1e6,
        depth_mmag=depth_mmag,
        t14_h=t14_h,
        t23_h=t23_h,
        ingress_min=ingress_m,
        snr_per_frame=snr_frame,
        noise_per_frame_mmag=sigma_flux * 1e3,
        depth_over_sigma=depth_mmag / (sigma_flux * 1e3),
        max_safe_exposure_s=calculate_max_exposure(star_mag),
        frames_in_transit=n_in,
        bins_in_transit=n_bin_in,
        binned_sigma_mmag=binned_sigma_mmag,
        grazing=b_imp >= 1 - k,
    )
