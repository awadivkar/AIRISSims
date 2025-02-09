import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import filedialog, messagebox
from astropy.io import fits
from scipy.ndimage import gaussian_filter

# Default Constants
DEFAULTS = {
    "Pixel Size (um)": 3.76e-6 * 3,
    "Sensor Width (px)": int(9568 / 3),
    "Sensor Height (px)": int(6380 / 3),
    "Aperture": 0.111,
    "QE": 0.6,
    "Wavelength (nm)": 640,
    "Dark Current (e-)": 0.56,
    "Saturation Capacity (e-)": 51000,
    "Readout Noise (e-)": 1,
    "Field of View (deg)": 10,
    "Max Magnitude": 20,
    "Min Magnitude": 0,
    "Zero Point": 10,
    "PSF (sigma)": 3,
    "F_Vega": 2.55e-12,
    "Exposure Time": 10,
    "Num of Stars": 1000
}

def generate_image(params):
    image = np.zeros((int(params["Sensor Height (px)"]), int(params["Sensor Width (px)"])))
    x_positions = np.random.uniform(0, params["Sensor Width (px)"], int(params["Num of Stars"]))
    y_positions = np.random.uniform(0, params["Sensor Height (px)"], int(params["Num of Stars"]))
    magnitudes = np.random.uniform(params["Min Magnitude"], params["Max Magnitude"], int(params["Num of Stars"]))
    
    for x, y, mag in zip(x_positions, y_positions, magnitudes):
        flux = 10**(-0.4 * (mag - params["Zero Point"]))
        photons = flux * params["Exposure Time"]
        electrons = np.random.poisson(photons * params["QE"])
        add_psf(image, x, y, electrons, params["PSF (sigma)"])
    
    dark_noise = np.random.poisson(params["Dark Current (e-)"] * params["Exposure Time"], image.shape)
    readout_noise = np.random.normal(params["Readout Noise (e-)"], 1.5, image.shape).astype(int)
    image += dark_noise + readout_noise
    image = np.clip(image, 0, params["Saturation Capacity (e-)"]).astype(int)
    return image

def add_psf(image, x, y, flux, sigma):
    size = int(6 * sigma)
    y_indices, x_indices = np.meshgrid(np.arange(-size//2, size//2+1), np.arange(-size//2, size//2+1), indexing='ij')
    psf = np.exp(-(x_indices**2 + y_indices**2) / (2 * sigma**2))
    psf /= psf.sum()
    ix, iy = int(y), int(x)
    if 0 <= ix < image.shape[0] and 0 <= iy < image.shape[1]:
        x_start, x_end = max(0, ix - size//2), min(image.shape[0], ix + size//2+1)
        y_start, y_end = max(0, iy - size//2), min(image.shape[1], iy + size//2+1)
        sub_psf = psf[:x_end-x_start, :y_end-y_start]
        image[x_start:x_end, y_start:y_end] += flux * sub_psf

def save_image(image, filename, format):
    if format == 'png':
        plt.imsave(filename, image, cmap='gray')
    elif format == 'fits':
        fits.writeto(filename, image, overwrite=True)

def run_simulation():
    params = {key: float(entries[key].get()) if entries[key].get() else DEFAULTS[key] for key in DEFAULTS}
    image = generate_image(params)
    plt.imshow(image, cmap='gray', origin='lower')
    plt.colorbar(label='Electron Count')
    plt.title('Simulated CMOS Image')
    plt.show()

def save_file():
    filename = filedialog.asksaveasfilename(defaultextension=".png", filetypes=[("PNG Image", "*.png"), ("FITS File", "*.fits")])
    if filename:
        format = "fits" if filename.endswith(".fits") else "png"
        image = generate_image({key: float(entries[key].get()) if entries[key].get() else DEFAULTS[key] for key in DEFAULTS})
        save_image(image, filename, format)
        messagebox.showinfo("Success", f"Image saved as {filename}")

# GUI Setup
root = tk.Tk()
root.title("CMOS Image Simulation GUI")
entries = {}

for i, (key, value) in enumerate(DEFAULTS.items()):
    tk.Label(root, text=key).grid(row=i, column=0)
    entries[key] = tk.Entry(root)
    entries[key].grid(row=i, column=1)
    entries[key].insert(0, str(value))

tk.Button(root, text="Run Simulation", command=run_simulation).grid(row=len(DEFAULTS), column=0)
tk.Button(root, text="Save Image", command=save_file).grid(row=len(DEFAULTS), column=1)

root.mainloop()
