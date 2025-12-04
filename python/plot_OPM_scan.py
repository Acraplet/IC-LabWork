#!/usr/bin/env python3
"""
Read OPM_meas_XX_Y.txt files from../results/ and plot:

1) Mean power (nW) vs pulse width (ps) with sigma as error bar.
2) Mean power (nW) vs PD current (pA) with sigma as error bar.

Only include:
    measurement numbers 6–19 (inclusive) AND 21.

File format example:

File: data_powermeter_10_1.csv
Measurement number (XX): 10
Run number (Y): 1
Start of measurement: 2025-11-27 17:44:11
Sample interval [s]: 0.1
Number of samples: 244
Total duration [s]: 24.300000
Wavelength [nm]: 401.0
Pulse width [ps]: 850

PD current [pA]: 8.63

Power statistics (nW):
  Mean: 7.881553e+02
  Std:  1.058115e+00
  Min:  7.851700e+02
  Max:  7.909800e+02...
"""

from pathlib import Path
import re

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Configuration: which measurement numbers to use
# ---------------------------------------------------------------------
ALLOWED_MEAS_NUMBERS = range(760, 910, 10) #set(range(6, 20)) | {21}


def parse_result_file(path: Path):
    """
    Parse one OPM_meas_XX_Y.txt file.

    Returns:
        dict with keys:
            'meas_no', 'run_no', 'pulse_width_ps',
            'pd_current_pA', 'mean_power_nW', 'sigma_power_nW'
        or None if parsing fails.
    """

    text = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    meas_no = None
    run_no = None
    pulse_width_ps = None
    pd_current_pA = None
    mean_power_nW = None
    sigma_power_nW = None

    for line in text:
        s = line.strip()

        if s.startswith("Measurement number (XX):"):
            # Measurement number (XX): 10
            try:
                meas_no = int(s.split(":")[1])
            except Exception:
                pass

        elif s.startswith("Run number (Y):"):
            # Run number (Y): 1
            try:
                run_no = s.split(":")[1].strip()
            except Exception:
                pass

        elif s.startswith("Pulse width [ps]:"):
            # Pulse width [ps]: 850
            try:
                pulse_width_ps = float(s.split(":")[1])
            except Exception:
                pass

        elif s.startswith("PD current [pA]:"):
            # PD current [pA]: 8.63
            try:
                pd_current_pA = float(s.split(":")[1])
            except Exception:
                pass

        elif s.startswith("Mean:") and "Power statistics" in "\n".join(text):
            #   Mean: 7.881553e+02
            try:
                mean_power_nW = float(s.split(":")[1])
            except Exception:
                pass

        elif s.startswith("Std:"):
            #   Std:  1.058115e+00
            try:
                sigma_power_nW = float(s.split(":")[1])
            except Exception:
                pass

    # As a more robust backup for mean/std, search via regex:
    if mean_power_nW is None or sigma_power_nW is None:
        joined = "\n".join(text)
        m_mean = re.search(r"Mean:\s*([0-9.eE+-]+)", joined)
        m_std = re.search(r"Std:\s*([0-9.eE+-]+)", joined)
        if m_mean:
            mean_power_nW = float(m_mean.group(1))
        if m_std:
            sigma_power_nW = float(m_std.group(1))

    if None in (
        meas_no,
        run_no,
        pulse_width_ps,
        pd_current_pA,
        mean_power_nW,
        sigma_power_nW,
    ):
        print(f"Warning: could not fully parse {path.name}, skipping.")
        return None

    return {
        "meas_no": meas_no,
        "run_no": run_no,
        "pulse_width_ps": pulse_width_ps,
        "pd_current_pA": pd_current_pA,
        "mean_power_nW": mean_power_nW,
        "sigma_power_nW": sigma_power_nW,
    }


def main():
    results_dir = Path("../results/OPM_per_run")
    if not results_dir.is_dir():
        raise SystemExit(f"Directory not found: {results_dir}")

    # Collect data from all OPM_meas_*.txt files
    records = []
    for path in sorted(results_dir.glob("New_OPM_meas_*.txt")):
        r = parse_result_file(path)
        if r is None:
            continue

        if r["meas_no"] not in ALLOWED_MEAS_NUMBERS:
            # skip measurements not in 6–19 or 21
            continue

        records.append(r)

    if not records:
        raise SystemExit("No valid records found in../results with the given criteria.")

    # Convert to numpy arrays
    meas_no = np.array([r["meas_no"] for r in records])
    pulse_width = np.array([r["pulse_width_ps"] for r in records])
    pd_current = np.array([r["pd_current_pA"] for r in records])
    mean_power = np.array([r["mean_power_nW"] for r in records])
    sigma_power = np.array([r["sigma_power_nW"] for r in records])

    # Sort by pulse width for a nice plot
    order_pw = np.argsort(pulse_width)
    pulse_width_sorted = pulse_width[order_pw]
    mean_power_pw = mean_power[order_pw]
    sigma_power_pw = sigma_power[order_pw]

    # Sort by PD current
    order_pd = np.argsort(pd_current)
    pd_current_sorted = pd_current[order_pd]
    mean_power_pd = mean_power[order_pd]
    sigma_power_pd = sigma_power[order_pd]

    # ---------- NEW: write sorted by pulse width ----------
    out_path = Path("../results/new_pulse_power_pd_summary.txt")
    with out_path.open("w", encoding="utf-8") as f:
        f.write("pulse_width_ps mean_power_nW pd_current_pA\n")
        for pw, mp, pd in zip(pulse_width_sorted, mean_power_pw, pd_current):
            f.write(f"{pw:.0f} {mp:.6e} {pd:.6e}\n")
    # ------------------------------------------------------

    # ------------- Plot 1: mean power vs pulse width -----------------
    plt.figure(figsize=(7, 4))
    plt.errorbar(
        pulse_width_sorted,
        mean_power_pw,
        yerr=sigma_power_pw,
        fmt="o-",
        capsize=4,
        label="Mean power",
    )
    plt.xlabel("Pulse width [ps]")
    plt.ylabel("Mean power [nW]")
    plt.title("Mean power vs pulse width")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("../plots/summary/mean_power_vs_pulsewidth.pdf")
    plt.close()

    # ------------- Plot 2: mean power vs PD current ------------------
    plt.figure(figsize=(7, 4))
    plt.errorbar(
        pd_current_sorted,
        mean_power_pd,
        yerr=sigma_power_pd,
        fmt="o-",
        capsize=1,
        markersize = 1,
        label="Mean power",
    )
    plt.xlabel("PD current [pA]")
    plt.ylabel("Mean power [nW]")
    plt.title("Mean power vs PD current")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.yscale("log")
    plt.savefig("../plots/summary/mean_power_vs_PDcurrent.pdf")
    plt.savefig("../plots/summary/mean_power_vs_PDcurrent.png")
    plt.close()

    print("Created:")
    print("  ../plots/summary/mean_power_vs_pulsewidth.pdf")
    print("  ../plots/summary/mean_power_vs_PDcurrent.pdf")


if __name__ == "__main__":
    main()
