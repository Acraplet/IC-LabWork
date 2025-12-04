import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# --- Configuration ---
txt_file      = "../results/Monitor_PMT_volts_summary.txt"                     # Monitor PMT summary
txt_file_OPM  = "../results/new_pulse_power_pd_summary.txt"  # OPM summary
out_pdf       = "../plots/summary/MonitorPMT_mean_vs_OPM.pdf"
out3_pdf       = "../plots/summary/MonitorPMT_minmax_vs_pd_current.pdf"
out2_pdf       = "../plots/summary/MonitorPMT_mean_vs_pd_current.pdf"

# ------------------------------------------------------------------
# Read data
# ------------------------------------------------------------------
data      = np.loadtxt(txt_file,     skiprows=1)
data_OPM  = np.loadtxt(txt_file_OPM, skiprows=1)

# First dataset: Monitor PMT (from ROOT analysis)
pulse_width     = data[:, 0]
sum_vals        = data[:, 1] + 191045          # pedestal-corrected sum
max_vals        = data[:, 3] - data[:, 2]      # max - min
mean            = data[:, 4]                   # mean integrated value
std             = data[:, 5]                   # sigma on mean

# Second dataset: OPM + PD current
pulse_width_OPM = data_OPM[:, 0]
OPM_value       = data_OPM[:, 1]
pd_current      = data_OPM[:, 2]

# ------------------------------------------------------------------
# Build dataframes and match by pulse width
# ------------------------------------------------------------------
df1 = pd.DataFrame({
    "pulse_width": pulse_width,
    "sum_vals": sum_vals,
    "max_vals": max_vals,
    "mean": mean,
    "std": std,
})

df2 = pd.DataFrame({
    "pulse_width": pulse_width_OPM,
    "OPM_value": OPM_value,
    "pd_current": pd_current,
})

# Remove pulse_width == 0 (pedestal row)
df1 = df1[df1["pulse_width"] != 0]

# Inner join: keep only pulse widths present in both datasets
df = pd.merge(df1, df2, on="pulse_width", how="inner")

if df.empty:
    raise SystemExit("No overlapping pulse_width entries between Monitor PMT and OPM files.")

# ------------------------------------------------------------------
# Prepare arrays for plotting and fitting
# ------------------------------------------------------------------
x = df["OPM_value"].values             # OPM mean power (e.g. nW)
y = df["mean"].values                  # Monitor PMT mean integrated value
yerr = df["std"].values                # Monitor PMT sigma

# ------------------------------------------------------------------
# Fit a straight line: y = a * x + b
# ------------------------------------------------------------------
# Use weights = 1 / sigma^2 for a weighted fit, if std > 0
def fit_straight_line(x, y, yerr):
    weights = None
    valid_err = yerr > 0
    if np.all(valid_err):
        weights = 1.0 / (yerr ** 2)
        p, cov = np.polyfit(x, y, deg=1, w=weights, cov=True)
    else:
        # fall back to unweighted fit
        p, cov = np.polyfit(x, y, deg=1, cov=True)

    slope, intercept = p
    slope_err    = np.sqrt(cov[0, 0])
    intercept_err = np.sqrt(cov[1, 1])

    # Build line for plotting
    x_fit = np.linspace(np.min(x), np.max(x), 200)
    y_fit = slope * x_fit + intercept

    return x_fit, y_fit, slope, intercept,  slope_err, intercept_err

x_fit, y_fit, slope, intercept, slope_err, intercept_err = fit_straight_line(x, y, yerr)

print("Fit results (Monitor PMT mean vs OPM):")
print(f"  y = a * x + b")
print(f"  a = {slope:.4g} ± {slope_err:.2g}")
print(f"  b = {intercept:.4g} ± {intercept_err:.2g}")

# ------------------------------------------------------------------
# Plot: Monitor PMT mean vs OPM with error bars and linear fit
# ------------------------------------------------------------------
plt.figure(figsize=(7, 5))

plt.errorbar(
    x,
    y,
    yerr=yerr,
    fmt="o",
    capsize=3,
    markersize=4,
    label="Data (Monitor PMT)"
)

plt.plot(
    x_fit,
    y_fit,
    "k--",
    label=rf"Fit: $y = a x + b$" "\n"
          rf"$a = {slope:.3g} \pm {slope_err:.2g}$" "\n"
          rf"$b = {intercept:.3g} \pm {intercept_err:.2g}$"
)

plt.xlabel("OPM mean power [nW]", fontsize = 18)
plt.ylabel("Monitor PMT int. peak [V]", fontsize = 18)
plt.title("Monitor PMT linearity vs OPM", weight = "bold", fontsize = 20)
plt.grid(True, alpha=0.3)
plt.legend(fontsize = 16)
plt.tight_layout()

plt.savefig(out_pdf)
plt.close()
########################################
plt.figure(figsize=(7, 5))

plt.errorbar(
    df["pd_current"].values ,
    y,
    yerr=yerr,
    fmt="o",
    capsize=3,
    markersize=4,
    label="Data (Monitor PMT)"
)

x_fit, y_fit, slope, intercept, slope_err, intercept_err = fit_straight_line(df["pd_current"].values, y, yerr)

plt.plot(
    x_fit,
    y_fit,
    "k--",
    label=rf"Fit: $y = a x + b$" "\n"
          rf"$a = {slope:.3g} \pm {slope_err:.2g}$" "\n"
          rf"$b = {intercept:.3g} \pm {intercept_err:.2g}$"
)

plt.xlabel("PD current [pA]", fontsize = 18)
plt.ylabel("Monitor PMT int. peak [V]" , fontsize = 18)
plt.title("Monitor PMT linearity vs PD current", weight = "bold", fontsize = 20)
plt.grid(True, alpha=0.3)
plt.legend(fontsize = 16)
plt.tight_layout()

plt.savefig(out2_pdf)
plt.close()

plt.figure(figsize=(7, 5))

plt.errorbar(
    df["pd_current"].values ,
    df["max_vals"].values,
    yerr=yerr*0,
    fmt="o",
    capsize=3,
    markersize=4,
    label="Data (Monitor PMT)"
)

x_fit, y_fit, slope, intercept, slope_err, intercept_err = fit_straight_line(df["pd_current"].values, df["max_vals"].values, yerr*0)

plt.plot(
    x_fit,
    y_fit,
    "k--",
    label=rf"Fit: $y = a x + b$" "\n"
          rf"$a = {slope:.3g} \pm {slope_err:.2g}$" "\n"
          rf"$b = {intercept:.3g} \pm {intercept_err:.2g}$"
)

plt.xlabel("PD current [pA]", fontsize = 18)
plt.ylabel("Monitor PMT Max-Min [V]" , fontsize = 18)
plt.title("Monitor PMT linearity vs PD current", weight = "bold", fontsize = 20)
plt.grid(True, alpha=0.3)
plt.legend(fontsize = 16)
plt.tight_layout()

plt.savefig(out3_pdf)
plt.close()

print(f"Saved plot to {out_pdf}")
print(f"Saved plot to {out2_pdf}")
print(f"Saved plot to {out3_pdf}")
