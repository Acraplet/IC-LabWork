#!/usr/bin/env python3
"""
Process Thorlabs powermeter CSV files (Optical Parameter Monitor).

Usage:
    python analyse_powermeter.py data_powermeter_07_0.csv --pulse-width 900

written by chat gpt, corrected by acraplet
"""

import re
import argparse
from pathlib import Path
from datetime import datetime, timedelta

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.optimize import curve_fit

#------------- Run information ------------
pulse_width_dict = {
    1:  900,
    2:  900,
    3:  900,
    4:  900,
    5:  900,
    6:  890,
    7:  880,
    8:  870,
    9:  860,
    10: 850,
    11: 840,
    12: 830,
    13: 820,
    14: 810,
    15: 800,
    16: 790,
    17: 780,
    18: 770,
    19: 760,
    20: 900,
    21: 900,


    }

pd_current_dict = {
    900: 9.171,
    890: 9.048,
    880: 8.926,
    870: 8.926,
    860: 8.806,
    850: 8.630,
    840: 8.400,
    830: 8.344,
    820: 7.958,
    810: 7.591,
    800: 7.489,
    790: 7.591,
    780: 7.489,
    770: 7.389,
    760: 7.290,
    }

# attenuation_dict = {
#     }

# ----------- Gaussian model --------------

def gaussian(x, A, mu, sigma):
    return A * np.exp(-0.5 * ((x - mu) / sigma) ** 2)


# ----------- Filename parsing -------------

def parse_filename_info(csv_path: Path):
    """Extract measurement number XX and run number Y from data_powermeter_XX_Y.csv"""
    m = re.search(r"data_powermeter_and_monitor_(\d+)_([0-9A-Za-z]+)", csv_path.stem)
    if not m:
        raise ValueError(f"Cannot parse measurement/run numbers from file name: {csv_path.name}")
    meas_no = int(m.group(1))
    run_no = m.group(2)
    return meas_no, run_no


# ----------- CSV reader (comma-delimited) -------------

def read_powermeter_csv(csv_path: Path):
    """
    Read powermeter CSV file with comma delimiter, header like:

    Start of Measurement , 11/27/2025 17:30:38
    Sample Interval , 0.1 s...
    Wavelength / Responsivity , 401 nm

    Samples,Date (MM/dd/yyyy),Time of day (hh:mm:ss),Power (W)
    0, 11/27/2025, 17:30:38.266,0.0000010757
    1, 11/27/2025, 17:30:38.653,0.0000010753...
    """

    lines = csv_path.read_text(encoding="utf-8", errors="ignore").splitlines()

    start_datetime = None
    sample_interval_s = None
    wavelength_nm = None
    header_end_idx = None

    for i, line in enumerate(lines):
        s = line.strip()

        if s.startswith("Start of Measurement"):
            # Example: Start of Measurement , 11/27/2025 17:30:38
            parts = [p.strip() for p in line.split(",") if p.strip()]
            # last two tokens should be date and time
            dt_str = parts[-1]  # "11/27/2025 17:30:38"
            start_datetime = datetime.strptime(dt_str, "%m/%d/%Y %H:%M:%S")

        elif s.startswith("Sample Interval"):
            # Example: Sample Interval , 0.1 s
            parts = [p.strip() for p in line.split(",") if p.strip()]
            # one of the parts will be something like "0.1 s"
            for p in parts:
                # take the first float-like token
                token = p.split()[0]
                try:
                    sample_interval_s = float(token)
                    break
                except ValueError:
                    continue

        elif s.startswith("Wavelength") or s.startswith("Wavelength / Responsivity"):
            # Example: Wavelength / Responsivity , 401 nm
            parts = [p.strip() for p in line.split(",") if p.strip()]
            for p in parts:
                token = p.split()[0]
                try:
                    wavelength_nm = float(token)
                    break
                except ValueError:
                    continue

        elif s.startswith("Samples") and "Power (W)" in s:
            header_end_idx = i
            break

    if header_end_idx is None:
        raise RuntimeError("Could not find data header row (starts with 'Samples').")

    if start_datetime is None:
        raise RuntimeError("Start of Measurement not found / not parsed.")
    if sample_interval_s is None:
        raise RuntimeError("Sample Interval not found / not parsed.")
    if wavelength_nm is None:
        raise RuntimeError("Wavelength not found / not parsed.")

    # ---- numeric data ----
    data_lines = [l for l in lines[header_end_idx + 1:] if l.strip()]
    if not data_lines:
        raise RuntimeError("No data lines found after 'Samples' header.")

    # Read as strings first to be robust
    arr = np.genfromtxt(
        data_lines,
        delimiter=",",
        dtype=str,
    )

    if arr.ndim == 1:
        arr = arr.reshape(1, -1)

    if arr.shape[1] < 4:
        raise RuntimeError(f"Expected at least 4 columns in data, got {arr.shape[1]}.")

    sample_index_str = arr[:, 0]
    power_str = arr[:, -1]

    sample_index = []
    power_W = []

    for s_idx, s_pow in zip(sample_index_str, power_str):
        try:
            i = int(s_idx.strip())
            p = float(s_pow.strip())
        except ValueError:
            # skip malformed rows (e.g. blank lines or trailing text)
            continue
        sample_index.append(i)
        power_W.append(p)

    sample_index = np.array(sample_index, dtype=int)
    power_W = np.array(power_W, dtype=float)

    if sample_index.size == 0:
        raise RuntimeError("No valid numeric data rows found.")

    # time in seconds since start
    time_s = sample_index * sample_interval_s
    power_nW = power_W * 1e9

    datetimes = np.array(
        [start_datetime + timedelta(seconds=float(t)) for t in time_s]
    )

    meta = {
        "start_datetime": start_datetime,
        "sample_interval_s": float(sample_interval_s),
        "wavelength_nm": float(wavelength_nm),
        "total_duration_s": float(time_s[-1] - time_s[0]),
        "n_samples": len(sample_index),
    }
    data = {
        "sample_index": sample_index,
        "datetime": datetimes,
        "time_s": time_s,
        "power_W": power_W,
        "power_nW": power_nW,
    }
    return meta, data


# ----------- Analysis & plotting -------------

def analyse_file(csv_path: Path, pulse_width_ps: float = 900.0):



    meas_no, run_no = parse_filename_info(csv_path)

    try:
        pulse_width_ps = pulse_width_dict[meas_no]
    except:
        if (meas_no > 600):
            pulse_width_ps = meas_no
        else:
            print("There was no pre-written pulse width value, we assume 900ps")
            pulse_width_ps = 900

    pd_current_pA = pd_current_dict[pulse_width_ps]

    # attenuation = 0.0 if (meas_no < 35 or meas_no > XX)  else attenuation_dict[meas_no]

    meta, data = read_powermeter_csv(csv_path)

    powers = data["power_nW"]

    # Histogram & Gaussian fit
    n_bins = max(20, int(np.sqrt(len(powers))))
    hist_counts, bin_edges = np.histogram(powers, bins=n_bins, density=False)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])

    A0 = hist_counts.max()
    mu0 = powers.mean()
    sigma0 = powers.std(ddof=1) if len(powers) > 1 else 1.0

    try:
        popt, _ = curve_fit(gaussian, bin_centers, hist_counts,
                            p0=[A0, mu0, sigma0])
        A_fit, mu_fit, sigma_fit = popt
    except RuntimeError:
        A_fit, mu_fit, sigma_fit = np.nan, np.nan, np.nan

    start_str = meta["start_datetime"].strftime("%Y-%m-%d %H:%M:%S")
    title_common = (
        f"Run {run_no} (meas. {meas_no}) - λ={meta['wavelength_nm']:.1f} nm, "
        f"pulse width={pulse_width_ps:.0f} ps, PD current = {pd_current_pA} pA \n"
        f"Start: {start_str}, Δt={meta['sample_interval_s']} s, "
        f"N={meta['n_samples']}, T={meta['total_duration_s']:.1f} s"
    )

    # ---- text output ----
    txt_path = f"../results/New_OPM_meas_{meas_no}_{run_no}.txt"
    with open(txt_path, "w") as f:
        f.write(f"File: {csv_path.name}\n")
        f.write(f"Measurement number (XX): {meas_no}\n")
        f.write(f"Run number (Y): {run_no}\n")
        f.write(f"Start of measurement: {start_str}\n")
        f.write(f"Sample interval [s]: {meta['sample_interval_s']}\n")
        f.write(f"Number of samples: {meta['n_samples']}\n")
        f.write(f"Total duration [s]: {meta['total_duration_s']:.6f}\n")
        f.write(f"Wavelength [nm]: {meta['wavelength_nm']}\n")
        f.write(f"Pulse width [ps]: {pulse_width_ps}\n\n")
        f.write(f"PD current [pA]: {pd_current_pA}\n\n")

        f.write("Power statistics (nW):\n")
        f.write(f"  Mean: {powers.mean():.6e}\n")
        f.write(f"  Std:  {powers.std(ddof=1):.6e}\n")
        f.write(f"  Min:  {powers.min():.6e}\n")
        f.write(f"  Max:  {powers.max():.6e}\n\n")

        f.write("Gaussian fit to histogram of power (nW):\n")
        f.write("  Model: A * exp(-0.5 * ((x - mu)/sigma)^2)\n")
        f.write(f"  A:     {A_fit:.6e}\n")
        f.write(f"  mu:    {mu_fit:.6e}  [nW]\n")
        f.write(f"  sigma: {sigma_fit:.6e}  [nW]\n")

    # ---- plots & PDF ----
    pdf_path = f"../plots/per_run/New_OPM_meas_{meas_no}_{run_no}.pdf"
    with PdfPages(pdf_path) as pdf:
        # Power vs time
        fig1, ax1 = plt.subplots(figsize=(9, 6))
        ax1.plot(data["time_s"], powers, lw=1)
        ax1.set_xlabel("Time since start [s]")
        ax1.set_ylabel("Power [nW]")
        ax1.set_title("Power vs. time\n" + title_common)
        ax1.grid(True, alpha=0.3)
        pdf.savefig(fig1)
        plt.close(fig1)

        # Histogram
        fig2, ax2 = plt.subplots(figsize=(9,6))
        ax2.hist(powers, bins=n_bins, alpha=0.6, color="C0",
                 label="Data", density=False)
        x_fit = np.linspace(powers.min(), powers.max(), 400)
        if not np.isnan(mu_fit):
            y_fit = gaussian(x_fit, A_fit, mu_fit, sigma_fit)
            ax2.plot(x_fit, y_fit, "r-", lw=2, label="Gaussian fit")
            leg_label = (
                "Fit parameters:\n"
                f"μ = {mu_fit:.3e} nW\n"
                f"σ = {sigma_fit:.3e} nW\n"
                f"σ/μ = {sigma_fit/mu_fit * 100:.3f} %"
            )
            ax2.legend(title=leg_label, fontsize=9)
        else:
            ax2.legend()

        ax2.set_xlabel("Power [nW]")
        ax2.set_ylabel("Counts")
        ax2.set_title("Histogram of power\n" + title_common)
        ax2.grid(True, alpha=0.3)
        pdf.savefig(fig2)
        plt.close(fig2)

    print(f"Processed {csv_path.name}")
    print(f"  -> Plots saved to: {pdf_path}")
    print(f"  -> Info saved to:  {txt_path}")


# ----------- CLI -------------

def main():
    parser = argparse.ArgumentParser(
        description="Analyse Thorlabs powermeter CSV data files."
    )
    parser.add_argument("csv_file", type=str, help="Input CSV file")
    parser.add_argument(
        "--pulse-width",
        type=float,
        default=900.0,
        help="Laser pulse width in ps (default: 900 ps)",
    )
    args = parser.parse_args()

    analyse_file(Path(args.csv_file), pulse_width_ps=args.pulse_width)


if __name__ == "__main__":
    main()
