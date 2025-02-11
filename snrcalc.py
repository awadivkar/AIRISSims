import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox
from astropy.io import fits
from scipy.ndimage import gaussian_filter

# Default Constants (including new parameters)
params = {
    "Pixel Size (um)": 3.76e-6,
    "Sensor Width (px)": 9568,
    "Sensor Height (px)": 6380,
    "Aperture": 0.111,
    "QE": 0.6,
    "Wavelength (nm)": 640,
    "Dark Current (e-)": 0.56,
    "Saturation Capacity (e-)": 51000,
    "Readout Noise (e-)": 1,
    "Field of View (deg)": 10,
    "Max Magnitude": 20,
    "Min Magnitude": 12,
    "Zero Point": 18,
    "PSF (sigma)": 3,
    "Exposure Time": 10,
    "Num of Stars": 1000,
    # New parameters for moving exposures:
    "Trail Length (px)": 10,
    "Drift Angle (deg)": 0,
    # New parameters for cosmic rays:
    "Cosmic Ray Count": 5,
    "Cosmic Ray Max Length": 20,
    "Cosmic Ray Intensity (e-)": 5000,
    # New parameter for sky background:
    "Sky Background Rate (e-/px/s)": 0.1
}

sky_background_var = True

signal = params["Exposure Time"] * params["QE"] * 10**(-0.4 * (params["Min Magnitude"] - params["Zero Point"]))
noise = np.sqrt(signal + int(sky_background_var) * params["Sky Background Rate (e-/px/s)"] + params["Dark Current (e-)"] * params["Exposure Time"] + params["Readout Noise (e-)"]**2)
snr = signal / noise
print("Expected Signal: ", signal, "Expected Noise: ", noise, "Expected SNR: ", snr)

x = np.linspace(0, 20, 100)
y = params["Exposure Time"] * params["QE"] * 10**(-0.4 * (x - params["Zero Point"]))
y = y / noise
plt.plot(x, y)
plt.yscale("log")
plt.axhline(y=3, color='r', linestyle='--')
plt.xlabel("Magnitude")
plt.ylabel("SNR")
plt.title("SNR vs. Magnitude")
plt.grid(True)
plt.show()