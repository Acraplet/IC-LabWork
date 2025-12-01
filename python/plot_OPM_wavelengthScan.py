#!/usr/bin/env python3
"""
Read OPM_meas_XX_Y.txt files from../results/ and plot:

Mean power (nW) vs wavelength (nm) with sigma as error bar
for measurement numbers 22–41 inclusive.

Each file is expected to contain lines like:

Wavelength [nm]: 401.0...
Gaussian fit to histogram of power (nW):
  mu:    7.880856e+02  [nW]
  sigma: 1.065547e+00  [nW]
"""

from pathlib import Path
import re

import numpy as np
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------
# Configuration: measurement numbers to keep
# ---------------------------------------------------------------------
ALLOWED_MEAS_NUMBERS = set(range(22, 42))  # 22–41 inclusive


def parse_result_file(path: Path):
    """
    Parse one OPM_meas_XX_Y.txt file.

    Returns:
        dict with keys:
            'meas_no', 'run_no',
            'wavelength_nm', 'mean_power_nW', 'sigma_power_nW'
        or None on failure.
    """
    lines = path.read_text(encoding="utf-8", errors="ignore").splitlines()

    meas_no = None
    run_no = None
    wavelength_nm = None

    # We'll use the Gaussian fit parameters (mu, sigma) as "power" and its error.
    gauss_mu = None
    gauss_sigma = None

    for line in lines:
        s = line.strip()

        if s.startswith("Measurement number (XX):"):
            try:
                meas_no = int(s.split(":")[1])
            except Exception:
                pass

        elif s.startswith("Run number (Y):"):
            try:
                run_no = s.split(":")[1].strip()
            except Exception:
                pass

        elif s.startswith("Wavelength [nm]:"):
            try:
                wavelength_nm = float(s.split(":")[1])
            except Exception:
                pass

        elif s.startswith("mu:") and "[nW]" in s:
            # Example: "  mu:    7.880856e+02  [nW]"
            try:
                # take the number before "[nW]"
                before_units = s.split("[nW]")[0]
                gauss_mu = float(before_units.split(":")[1])
            except Exception:
                pass

        elif s.startswith("sigma:") and "[nW]" in s:
            # Example: "  sigma: 1.065547e+00  [nW]"
            try:
                before_units = s.split("[nW]")[0]
                gauss_sigma = float(before_units.split(":")[1])
            except Exception:
                pass

    # As a backup, also try regex for mu/sigma in case spacing changes:
    if gauss_mu is None or gauss_sigma is None:
        joined = "\n".join(lines)
        m_mu = re.search(r"mu:\s*([0-9.eE+-]+)\s*\[nW\]", joined)
        m_sigma = re.search(r"sigma:\s*([0-9.eE+-]+)\s*\[nW\]", joined)
        if m_mu:
            gauss_mu = float(m_mu.group(1))
        if m_sigma:
            gauss_sigma = float(m_sigma.group(1))

    if None in (meas_no, run_no, wavelength_nm, gauss_mu, gauss_sigma):
        print(f"Warning: could not fully parse {path.name}, skipping.")
        return None

    return {
        "meas_no": meas_no,
        "run_no": run_no,
        "wavelength_nm": wavelength_nm,
        "mean_power_nW": gauss_mu,
        "sigma_power_nW": gauss_sigma,
    }


def main():
    results_dir = Path("../results")
    if not results_dir.is_dir():
        raise SystemExit(f"Directory not found: {results_dir}")

    records = []
    for path in sorted(results_dir.glob("OPM_meas_*.txt")):
        r = parse_result_file(path)
        if r is None:
            continue

        if r["meas_no"] not in ALLOWED_MEAS_NUMBERS:
            continue

        records.append(r)

    if not records:
        raise SystemExit("No valid records found in../results for meas 22–41.")

    # Convert to numpy arrays
    wavelengths = np.array([r["wavelength_nm"] for r in records])
    mean_power = np.array([r["mean_power_nW"] for r in records])
    sigma_power = np.array([r["sigma_power_nW"] for r in records])

    # Sort by wavelength for a nicer plot
    order = np.argsort(wavelengths)
    wavelengths = wavelengths[order]
    mean_power = mean_power[order]
    sigma_power = sigma_power[order]

    # -------- Plot: mean power vs wavelength --------
    plt.figure(figsize=(7, 4))
    plt.errorbar(
        wavelengths,
        mean_power,
        yerr=sigma_power,
        fmt="o-",
        capsize=4,
        label="Gaussian mean power",
    )
    plt.xlabel("Wavelength [nm]")
    plt.ylabel("Mean power [nW]")
    plt.title("Mean power vs wavelength\nPulse width = 900 ps, PD current = 9.171 pA")
    plt.grid(True, alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig("../plots/summary/power_vs_wavelength.pdf")
    plt.close()

    print("Created: ../plots/summary/power_vs_wavelength.pdf")


if __name__ == "__main__":
    main()
