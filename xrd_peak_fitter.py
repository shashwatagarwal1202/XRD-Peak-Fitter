# -*- coding: utf-8 -*-
"""
Interactive XRD Fitting
Voigt peaks + automatic residual-shoulder model
(Voigt shoulders with strict constraints + MANUAL PEAK FITTING RETAINED
 + HARD VOIGT TAIL TRUNCATION + SPLINE BACKGROUND + MID-SCALE POLY PROTECTION)
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.optimize import curve_fit
from scipy.signal import find_peaks, savgol_filter
from scipy.special import wofz
from scipy.interpolate import UnivariateSpline
from tkinter import Tk, filedialog
import os
import re
__author__ = "Shashwat Agarwal"
__version__ = "1.0.0"
# SPECTRUM TYPE = XRD


# ============================================================
# MODELS
# ============================================================

PEAK_PROMINENCE_FACTOR = 0.30
X_LABEL = "2θ"
Y_LABEL = "Intensity"
TRUNC_K = 5.0   # hard cutoff in units of FWHM
DISABLE_SHOULDERS = True
USER_EXCLUDE_X = []

def voigt(x, A, x0, sigma, gamma):
    z = ((x - x0) + 1j * gamma) / (sigma * np.sqrt(2))
    return A * np.real(wofz(z)) / (sigma * np.sqrt(2 * np.pi))

def voigt_truncated(x, A, x0, sigma, gamma):
    fwhm = 0.5346*(2*gamma) + np.sqrt(
        0.2166*(2*gamma)**2 + (2.3548*sigma)**2
    )
    mask = np.abs(x - x0) <= TRUNC_K * fwhm
    y = np.zeros_like(x)
    if np.any(mask):
        y[mask] = voigt(x[mask], A, x0, sigma, gamma)
    return y


def multi_peak_model(x, params):

    y = np.zeros_like(x)

    for A, x0, sigma, gamma in params:

        y += voigt_truncated(
            x,
            A,
            x0,
            sigma,
            gamma
        )

    return y


def peaks_only_model(x, *p):
    peaks = np.array(p).reshape(-1, 4)
    return multi_peak_model(x, peaks)


# ============================================================
# SPLINE BACKGROUND
# ============================================================

def spline_background(x, y, peak_mask, s_factor=0.8):
    x_bg = x[~peak_mask]
    y_bg = y[~peak_mask]
    s = s_factor * len(x_bg)
    spline = UnivariateSpline(x_bg, y_bg, s=s)
    return spline(x)


# ============================================================
# MID-SCALE POLYNOMIAL BACKGROUND (PROTECTION TERM)
# ============================================================

def poly_background(x,y, peak_mask, degree=3):
    x_bg = x[~peak_mask]
    y_bg = y[~peak_mask]
    coeffs = np.polyfit(x_bg,y_bg, degree)
    return np.polyval(coeffs, x)

# ============================================================
# FILE HANDLING
# ============================================================

def pick_files():
    root = Tk()
    root.withdraw()
    paths = filedialog.askopenfilenames(
        title="Select XRD files",
        filetypes=[("Data files", "*.txt *.dat *.csv *.asc *.chi"),
                   ("All files", "*.*")]
    )
    root.destroy()
    if not paths:
        raise RuntimeError("No files selected.")
    return list(paths)

def ask_save_location():

    root = Tk()
    root.withdraw()

    save_path = filedialog.asksaveasfilename(
        title="Save Results",
        defaultextension=".csv",
        filetypes=[
            ("CSV files", "*.csv"),
            ("Excel files", "*.xlsx"),
            ("NPZ files", "*.npz"),
            ("JSON files", "*.json"),
            ("Text files", "*.txt"),
            ("All files", "*.*")
        ]
    )

    root.destroy()

    return save_path

def load_any(path):
    rows = []
    with open(path, "r", errors="ignore") as f:
        for line in f:
            line = line.replace(",", " ").replace("\t", " ")
            nums = []
            for t in line.split():
                try:
                    nums.append(float(t))
                except:
                    pass
            if len(nums) >= 2:
                rows.append(nums[:2])
    data = np.array(rows)
    return data[:,0], data[:,1]

# ============================================================
# SAFE PEAK FIT
# ============================================================

def fit_peaks_safely(x, y_corr, window, peaks_guess):
    p0 = peaks_guess.flatten()
    lower, upper = [], []

    for A, x0, s, g in peaks_guess:
        lower += [0.0, x.min(), 1e-6, 1e-6]

        upper += [
            2 * A if A > 0 else np.inf,
            x.max(),
            window,
            window
        ]

    popt, _ = curve_fit(
        peaks_only_model,
        x,
        y_corr,
        p0=p0,
        bounds=(lower, upper),
        maxfev=15000
    )

    return np.array(popt).reshape(-1,4)

def exclude_peaks(
    x,
    y,
    peaks_fit,
    bg_est,
    poly_bg
):
    print("Click near the TOP of peaks to exclude. Enter/Right-click to finish.")

    fig, ax = plt.subplots(figsize=(9, 5))

    model = peaks_only_model(x, *peaks_fit.flatten()) + bg_est + poly_bg

    ax.scatter(x, y, s=12)
    ax.plot(x, model, lw=2)

    # Mark peak centers
    for i, (_, x0, _, _) in enumerate(peaks_fit):
        ax.axvline(x0, color='gray', alpha=0.3)
        ax.text(x0, max(y)*0.9, str(i), rotation=90,
                ha='center', va='top', fontsize=9)

    plt.title("Click near peak TOP to remove (numbers shown)")

    clicked = plt.ginput(n=-1, timeout=0)
    global USER_EXCLUDE_X
    USER_EXCLUDE_X.extend([c[0] for c in clicked])

    exclude_peaks.last_clicks = clicked
    global USER_EXCLUDE_MASK
    plt.close(fig)

    if not clicked:
        return peaks_fit

    exclude = set()

    # tolerances
    x_tol = 0.02 * (x.max() - x.min())
    y_tol = 0.15 * np.max(y)

    for xc, yc in clicked:

        dx = np.abs(peaks_fit[:,1] - xc)

        idx = np.argmin(dx)

        # Check if click is close enough in X
        if dx[idx] > x_tol:
            continue

        # Check if click is near peak height
        peak_y = np.interp(peaks_fit[idx,1], x, model)

        if abs(peak_y - yc) > y_tol:
            continue

        exclude.add(idx)

    if not exclude:
        print("No valid peaks selected.")
        return peaks_fit

    exclude = sorted(exclude)
    print("Removed peaks:", exclude)

    return np.delete(peaks_fit, exclude, axis=0)

def apply_user_exclusion(x, y, width_frac=0.01):

    if not USER_EXCLUDE_X:
        return x, y

    mask = np.ones_like(x, dtype=bool)
    x_range = x.max() - x.min()
    half_w = width_frac * x_range

    for xc in USER_EXCLUDE_X:
        mask &= np.abs(x - xc) > half_w

    return x[mask], y[mask]


def build_exclusion_mask(x, clicks, width_frac=0.01):
    """
    Build mask for excluded x-regions around user clicks
    width_frac = fraction of x-range to exclude per click
    """
    mask = np.ones_like(x, dtype=bool)

    if not clicks:
        return mask

    x_range = x.max() - x.min()
    half_width = width_frac * x_range

    for xc, _ in clicks:
        mask &= np.abs(x - xc) > half_width

    return mask

exclude_peaks.last_clicks = []

def initial_exclusion(x, y):

    print("Initial cleanup: Click anywhere to EXCLUDE regions. Enter to continue.")

    fig, ax = plt.subplots(figsize=(9,5))
    ax.scatter(x, y, s=10)
    plt.title("Initial exclusion (click junk/artifacts)")
    clicks = plt.ginput(n=-1, timeout=0)
    plt.close(fig)

    if not clicks:
        return x, y

    mask = np.ones_like(x, dtype=bool)

    x_range = x.max() - x.min()
    half_w = 0.01 * x_range   # exclusion width

    for xc, _ in clicks:
        mask &= np.abs(x - xc) > half_w

    return x[mask], y[mask]

# ============================================================
# MAIN
# ============================================================

def run_analysis():
    paths = pick_files()

    all_results = []
    all_dfs = []

    for path in paths:

        USER_EXCLUDE_X.clear()

        print(f"\nProcessing: {os.path.basename(path)}")

        x, y = load_any(path)
        idx = np.argsort(x)
        x, y = x[idx], y[idx]

        # ---- INITIAL USER CLEANUP (BEFORE EVERYTHING) ----
        x, y = initial_exclusion(x, y)

        # ---- Build user exclusion mask (ONCE) ----
        if USER_EXCLUDE_X:
            user_exclude_mask = build_exclusion_mask(x, [(v, 0) for v in USER_EXCLUDE_X])
        else:
            user_exclude_mask = np.ones_like(x, dtype=bool)
        # ---- Apply exclusion ----
        x = x[user_exclude_mask]
        y = y[user_exclude_mask]

        # ---- Apply user exclusion BEFORE background ----
        x, y = apply_user_exclusion(x, y)

        # ---- Rough peak mask ----
        y_s = savgol_filter(y, 31, 3)
        rough_peaks, _ = find_peaks(y_s, prominence=0.2 * np.std(y_s))

        peak_mask = np.zeros_like(y, dtype=bool)
        w = int(0.02 * len(y))
        for p in rough_peaks:
            peak_mask[max(0, p - w):min(len(y), p + w)] = True

        # ---- Apply user exclusion to mask ----
        if user_exclude_mask is not None:

            # ---- Spline background ----
            bg_est = spline_background(
                x,
                y,
                peak_mask
            )

            bg_est -= np.median(bg_est)

            y_corr = y - bg_est
            # ---- Apply user exclusion to data ----
            # if user_exclude_mask is not None:
            #     x = x[user_exclude_mask]
            #     y = y[user_exclude_mask]
            #     y_corr = y_corr[user_exclude_mask]

        # ---- Polynomial mid-scale background (NEW) ----
        poly_bg = poly_background(x, y_corr, peak_mask, degree=3)
        y_corr -= poly_bg

        # ---- Peak detection ----
        y_s = savgol_filter(y_corr, 31, 3)
        peaks, _ = find_peaks(
            y_s,
            prominence=PEAK_PROMINENCE_FACTOR * np.std(y_s)
        )

        window = (x.max() - x.min()) / (6 * max(len(peaks), 1))
        window = max(window, 1e-4 * (x.max() - x.min()))
        peak_params = []

        for p in peaks:
            A0 = max(
                y_corr[p],
                0.05 * np.std(y_corr)
            )

            peak_params.append(
                [
                    A0,
                    x[p],
                    window / 4,
                    window / 4
                ]
            )

        # ---- Initial peak fit ----
        peaks_fit = fit_peaks_safely(x, y_corr, window, np.array(peak_params))

        # ---- Manual peak addition ----
        peaks_fit = exclude_peaks(
            x,
            y_corr,
            peaks_fit,
            bg_est,
            poly_bg
        )

        # ---- Shoulders ----
        shoulder_fit = np.zeros_like(x)

        print("(Optional) Click to add peaks. Right click / Enter to finish.")

        fig, ax = plt.subplots(figsize=(9, 5))
        ax.scatter(x, y, s=12)
        ax.plot(x, peaks_only_model(x, *peaks_fit.flatten()) + bg_est + poly_bg)
        plt.title(os.path.basename(path))
        clicked = plt.ginput(n=-1, timeout=0)
        plt.close(fig)

        for (x_click, _) in clicked:
            idx_c = np.argmin(
                np.abs(x - x_click)
            )

            A = max(
                y_corr[idx_c],
                0.05 * np.std(y_corr)
            )

            sigma = window / 4

            peaks_fit = np.vstack(
                [
                    peaks_fit,
                    [A, x_click, sigma, sigma]
                ]
            )

        if len(clicked) > 0:
            peaks_fit = fit_peaks_safely(
                x,
                y_corr,
                window,
                peaks_fit
            )

        # ---- Shoulders ----
        if not DISABLE_SHOULDERS:

            y_voigt = peaks_only_model(x, *peaks_fit.flatten()) + bg_est + poly_bg
            residual = y - y_voigt
            shoulder_fit = np.zeros_like(x)

            mean_sigma = np.mean(peaks_fit[:, 2])
            mean_gamma = np.mean(peaks_fit[:, 3])

            r_s = savgol_filter(residual, 21, 3)
            thr = 0.05 * np.max(np.abs(residual))
            mask = np.abs(r_s) > thr

            regions, current = [], []
            for i, flag in enumerate(mask):
                if flag:
                    current.append(i)
                else:
                    if len(current) > 10:
                        regions.append(current)
                    current = []
            if len(current) > 10:
                regions.append(current)

            for reg in regions:
                xr = x[reg]
                rr = residual[reg]
                x0 = np.average(xr, weights=np.abs(rr))
                sigma_s = max((xr.max() - xr.min()) / 2.5, 2.0 * mean_sigma)
                gamma_s = max(mean_gamma, 0.5 * sigma_s)
                V = voigt_truncated(x, 1.0, x0, sigma_s, gamma_s)
                A = max(0.0, np.dot(V, residual) / np.dot(V, V))
                shoulder_fit += A * V

        if DISABLE_SHOULDERS:
            shoulder_fit = np.zeros_like(x)

        # ---- Final plot ----
        total_fit = peaks_only_model(x, *peaks_fit.flatten()) \
                    + bg_est + poly_bg + shoulder_fit

        residuals = y - total_fit
        percent_residual = 100 * residuals / np.maximum(y, 1e-6)

        print("\nFit Quality")
        print(
            f"Mean abs % error = "
            f"{np.mean(np.abs(percent_residual)):.3f}%"
        )
        rmse = np.sqrt(
            np.mean(residuals ** 2)
        )

        ss_res = np.sum(
            residuals ** 2
        )

        ss_tot = np.sum(
            (y - np.mean(y)) ** 2
        )

        r2 = 1 - ss_res / ss_tot

        print(f"RMSE = {rmse:.6f}")
        print(f"R²   = {r2:.8f}")
        max_res = np.max(np.abs(residuals))

        print(f"Max residual = {max_res:.4f}")
        x_dense = np.linspace(x.min(), x.max(), 4000)
        bg_dense = np.interp(x_dense, x, bg_est)
        poly_dense = np.interp(x_dense, x, poly_bg)
        y_voigt_dense = peaks_only_model(x_dense, *peaks_fit.flatten()) + bg_dense + poly_dense
        shoulder_dense = np.interp(x_dense, x, shoulder_fit)

        fig, (ax1, ax2) = plt.subplots(
            2,
            1,
            figsize=(10, 7),
            sharex=True,
            gridspec_kw={"height_ratios": [3, 1]}
        )

        # =========================
        # MAIN FIT PLOT
        # =========================

        ax1.scatter(
            x,
            y,
            s=12,
            color='black',
            label='Data'
        )

        ax1.plot(
            x_dense,
            y_voigt_dense + shoulder_dense,
            lw=2.5,
            color='red',
            label='Total Fit'
        )

        ax1.plot(
            x_dense,
            y_voigt_dense,
            '--',
            lw=2,
            color='green',
            label='Voigt Peaks'
        )

        ax1.plot(
            x_dense,
            shoulder_dense,
            ':',
            lw=2,
            color='purple',
            label='Shoulders'
        )

        ax1.plot(
            x_dense,
            bg_dense + poly_dense,
            '--',
            lw=2,
            color='blue',
            label='Background'
        )

        ax1.legend()

        ax1.set_title(
            f"Final Fit: {os.path.basename(path)}"
        )

        ax1.text(
            0.01,
            0.98,
            "Spectrum: XRD",
            transform=ax1.transAxes,
            verticalalignment="top"
        )

        # =========================
        # RESIDUAL PLOT
        # =========================

        ax2.axhline(
            0,
            color='black',
            lw=1
        )

        ax2.scatter(
            x,
            residuals,
            s=10,
            color='red'
        )

        ax2.set_ylabel("Residual")

        ax1.set_ylabel(Y_LABEL)

        ax2.set_xlabel(X_LABEL)

        fig.tight_layout()
        plt.show()

        # ---- Save parameters ----
        output_peaks = []

        for peak_num, (A, x0, s, g) in enumerate(
                peaks_fit,
                start=1
        ):
            fwhm = 0.5346 * (2 * g) + np.sqrt(
                0.2166 * (2 * g) ** 2 + (2.3548 * s) ** 2
            )

            output_peaks.append(
                [
                    peak_num,
                    x0,
                    abs(A),
                    abs(s),
                    abs(g),
                    fwhm
                ]
            )

        output_peaks = np.array(output_peaks)

        df = pd.DataFrame(
            output_peaks,
            columns=[
                "Peak Number",
                "Center",
                "Amplitude",
                "Sigma",
                "Gamma",
                "FWHM"
            ]
        )

        df["RMSE"] = rmse
        df["R2"] = r2
        df["Max Residual"] = max_res

        all_dfs.append(
            (os.path.basename(path), df.copy())
        )

        output_peaks = output_peaks[np.argsort(output_peaks[:, 1])]
        all_results.append((os.path.basename(path), output_peaks))

    save_path = ask_save_location()

    if not save_path:
        print("\nResults not saved.")
        exit()

    ext = os.path.splitext(save_path)[1].lower()

    if ext == "":
        ext = ".csv"
        save_path += ".csv"

    if ext == ".npz":

        npz_data = {}

        for fname, peaks in all_results:
            key = os.path.splitext(fname)[0]

            npz_data[key] = peaks.astype(np.float32)

        np.savez_compressed(
            save_path,
            **npz_data
        )

        print(
            f"\nSaved NPZ:\n{save_path}"
        )


    elif ext == ".csv":

        combined_df = pd.concat(
            [df for _, df in all_dfs],
            ignore_index=True
        )

        combined_df.to_csv(
            save_path,
            index=False
        )

        print(
            f"\nSaved CSV:\n{save_path}"
        )


    elif ext == ".xlsx":

        with pd.ExcelWriter(save_path) as writer:

            for fname, df in all_dfs:
                sheet_name = re.sub(
                    r'[:\\/*?\[\]]',
                    '_',
                    os.path.splitext(fname)[0]
                )[:31]

                df.to_excel(
                    writer,
                    sheet_name=sheet_name,
                    index=False
                )

        print(
            f"\nSaved Excel:\n{save_path}"
        )

    elif ext == ".json":

        combined_df = pd.concat(
            [df for _, df in all_dfs],
            ignore_index=True
        )

        combined_df.to_json(
            save_path,
            orient="records",
            indent=4
        )

        print(
            f"\nSaved JSON:\n{save_path}"
        )

    elif ext in [".txt", ".dat"]:

        combined_df = pd.concat(
            [df for _, df in all_dfs],
            ignore_index=True
        )

        combined_df.to_csv(
            save_path,
            sep="\t",
            index=False
        )

        print(
            f"\nSaved TXT/DAT:\n{save_path}"
        )

    else:

        raise RuntimeError(
            f"Unsupported format: {ext}"
        )

if __name__ == "__main__":
    run_analysis()