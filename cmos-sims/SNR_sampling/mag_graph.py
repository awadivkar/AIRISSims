import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
import pandas as pd
from scipy.optimize import curve_fit

matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = ["Inter", "DejaVu Sans", "Arial", "Liberation Sans"]
# matplotlib.rcParams["font.weight"] = "bold"

# Default Constants (including new parameters)
DEFAULTS = {
    "Pixel Size (um)": 3.76e-6,
    "Sensor Width (px)": 50,
    "Sensor Height (px)": 50,
    "Aperture": 0.111,
    "QE": 0.6,
    "Wavelength (nm)": 640,
    "Dark Current (e-)": 0.56,
    "Saturation Capacity (e-)": 51000,
    "Readout Noise (e-)": 1,
    "Width Field of View (deg)": 10,
    "Min Magnitude": 10,
    "Max Magnitude": 20,
    "Zero Point": 18,
    "PSF (sigma)": 3,
    "Exposure Time": 1,
}

OPTIONAL_PARAMS = {
    "Trail Length (px)": 10,
    "Drift Angle (deg)": 0,
    "Cosmic Ray Count": 5,
    "Cosmic Ray Max Length": 20,
    "Cosmic Ray Intensity (e-)": 5000,
    "Sky Background Rate (e-/px/s)": 0.1,
}

var = 'SNR'
datapath = 'cmos-sims/SNR_sampling/mag_stats.csv'
datapath3 = 'cmos-sims/SNR_sampling/mag_stats(3x3).csv'

df = pd.read_csv(datapath)
df3 = pd.read_csv(datapath3)
df = df[df['Magnitude'] > 10]
df3 = df3[df3['Magnitude'] > 10]

def fit_func(x, a, b, c):
    return a * np.exp(-b * x) + c

# Fit the data to the function
popt, pcov = curve_fit(fit_func, df['Magnitude'], df[var])
popt3, pcov3 = curve_fit(fit_func, df3['Magnitude'], df3[var])
# Generate x values for the fitted curve
x_fit = np.linspace(df['Magnitude'].min(), df['Magnitude'].max(), 100)
x_fit3 = np.linspace(df3['Magnitude'].min(), df3['Magnitude'].max(), 100)
# Calculate the fitted y values
y_fit = fit_func(x_fit, *popt)
y_fit3 = fit_func(x_fit3, *popt3)
# Plot the data and the fitted curve

print(f"mean difference; {(np.mean(df[var] - df3[var])) / np.mean(df[var]) * 100:.2f}%")

fig, ax = plt.subplots(figsize=(8, 6))
ax.scatter(df['Magnitude'], df[var], marker='o', alpha=0.1, color='k', label=f'Simulated {var} vs Magnitude')
ax.scatter(df3['Magnitude'], df3[var], marker='o', alpha=0.1, color='b', label=f'3x3 Binned {var} vs Magnitude')
ax.plot(x_fit, y_fit, color='orange', label=f'Fitted Curve: {var} = {popt[0]:.3G} $\\times \ \exp(-{popt[1]:.2f} \\times Magnitude) - {-popt[2]:.2f}$')
ax.plot(x_fit3, y_fit3, color='blue', label=f'Fitted Curve (3x3 Binned): {var} = {popt3[0]:.3G} $\\times \ \exp(-{popt3[1]:.2f} \\times Magnitude) - {-popt3[2]:.2f}$')
# plt.axhline(y=10, color='g', linestyle='--', label='SNR = 10 Threshold')
plt.axhline(y=5, color='g', linestyle='--', label='SNR = 5 Threshold')
plt.axhline(y=3, color='r', linestyle='--', label='SNR = 3 Threshold')

ax.set_xlabel('Magnitude', fontsize=14)
ax.set_xbound(10, 16)
ax.set_ylabel('SNR  $ (\mu / \sigma)$', fontsize=14)
ax.set_yscale('log')
ax.set_title(f'{var} vs Magnitude', fontsize=14)
# ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
# ax.yaxis.get_major_formatter().set_scientific(True)
ax.grid(True, which='both', linestyle='--', linewidth=0.5)
ax.legend(fontsize=12)
plt.show()



# plt.hist(df[' Magnitude'], bins=200, color='blue', alpha=0.7)
# plt.xlabel('Magnitude')
# plt.ylabel('Frequency')
# plt.title('Histogram of Magnitude')
# plt.grid(True)
# plt.show()



# sky_background_var = False

# signal = DEFAULTS["Exposure Time"] * DEFAULTS["QE"] * 10**(-0.4 * (DEFAULTS["Min Magnitude"] - DEFAULTS["Zero Point"]))
# noise = np.sqrt(signal 
#                 + int(sky_background_var) * OPTIONAL_PARAMS["Sky Background Rate (e-/px/s)"] 
#                 + DEFAULTS["Dark Current (e-)"] * DEFAULTS["Exposure Time"] 
#                 + DEFAULTS["Readout Noise (e-)"]**2)
# snr = signal / noise
# print("Expected Signal: ", signal, "Expected Noise: ", noise, "Expected SNR: ", snr)

# x = np.linspace(0, 20, 100)
# y = DEFAULTS["Exposure Time"] * DEFAULTS["QE"] * 10**(-0.4 * (x - DEFAULTS["Zero Point"]))
# y = y / noise
# plt.plot(x, y)
# plt.yscale("log")
# plt.axhline(y=3, color='r', linestyle='--')
# plt.xlabel("Magnitude")
# plt.ylabel("SNR")
# plt.title("SNR vs. Magnitude")
# plt.grid(True)
# plt.show()

