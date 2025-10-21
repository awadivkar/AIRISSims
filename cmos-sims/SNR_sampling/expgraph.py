import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
import pandas as pd
from scipy.optimize import curve_fit

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
    "Exposure Time": 10,
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
datapath = 'expstats.csv'
df = pd.read_csv(datapath)
# df = df[df[' Magnitude'] > 10]

def fit_func(x, a, b, c):
    return a * np.exp(-b * x) + c

# Fit the data to the function
popt, pcov = curve_fit(fit_func, df['Exposure Time'], df[' ' + var])
# Generate x values for the fitted curve
x_fit = np.linspace(df['Exposure Time'].min(), df['Exposure Time'].max(), 100)
# Calculate the fitted y values
y_fit = fit_func(x_fit, *popt)
# Plot the data and the fitted curve


plt.figure(figsize=(10, 6))
plt.scatter(df['Exposure Time'], df[' ' + var], marker='o', alpha=0.1, color='b', label=f'{var} vs Exposure Time')
# plt.plot(x_fit, y_fit, color='orange', label=f'Fitted Curve: {var} = {popt[0]:.2f} * exp(-{popt[1]:.2f} * Exposure Time) + {popt[2]:.2f}')
# plt.axhline(y=10, color='g', linestyle='--', label='SNR = 10 Threshold')
# plt.axhline(y=5, color='y', linestyle='--', label='SNR = 5 Threshold')
# plt.axhline(y=3, color='r', linestyle='--', label='SNR = 3 Threshold')
plt.xlabel('Exposure Time')
plt.ylabel(var)
# plt.xscale('log')
# plt.yscale('log')
plt.title(f'{var} vs Exposure Time')
plt.grid(True)
plt.legend()
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

