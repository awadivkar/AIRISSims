import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.widgets import Slider
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from astropy.io import fits
from skimage.transform import resize
from PIL import Image, ImageTk

# Default Constants
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
    "Star Magnitude": 10,
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

## Helper Functions

def apply_binning(image, bin_size=3):
    """Bin the image in non-overlapping blocks of bin_size x bin_size pixels."""
    h, w = image.shape
    # Trim the image so that dimensions are divisible by bin_size:
    h_trim = h - (h % bin_size)
    w_trim = w - (w % bin_size)
    trimmed = image[:h_trim, :w_trim]
    # Reshape and sum over the binning blocks.
    binned = trimmed.reshape(h_trim // bin_size, bin_size, w_trim // bin_size, bin_size).sum(axis=(1,3))
    return binned

def add_cosmic_rays(image, num_rays=5, max_length=20, intensity=5000):
    ## Add simulated cosmic ray events as short bright streaks."""
    for _ in range(num_rays):
        # Choose a random starting pixel.
        start_x = np.random.randint(0, image.shape[1])
        start_y = np.random.randint(0, image.shape[0])
        # Random length and direction.
        length = np.random.randint(1, max_length)
        angle = np.random.uniform(0, 2 * np.pi)
        for i in range(length):
            x = int(start_x + i * np.cos(angle))
            y = int(start_y + i * np.sin(angle))
            if 0 <= x < image.shape[1] and 0 <= y < image.shape[0]:
                image[y, x] += intensity
    return image

def add_sky_background(image, background_rate, exposure_time):
    """Add a uniform sky background (in electrons) across the image."""
    # In a more detailed model, you’d need the sun’s elevation, atmospheric conditions, etc.
    # Check skybackcalc.py for a (presumably) more accurate calculation.
    return image + background_rate * exposure_time

def calculate_image_flux(image_magnitude, sensor_params):
    """Calculate appropriate photon flux per pixel for the imported image."""
    total_flux = 10 ** (-0.4 * (image_magnitude - sensor_params["Zero Point"]))
    total_photons = total_flux * sensor_params["Exposure Time"]  # Adjust based on exposure time
    total_photons *= sensor_params["QE"]  # Adjust based on quantum efficiency
    # total_photons *= 0.1 # Adjust based on \delta 30 nm bandpass filter
    pixel_flux = total_photons / np.sum(imported_image)  # Normalize across image pixels
    return pixel_flux

def add_imported_image(image, scale_factor=1.0, image_magnitude=10.0, sensor_params=None):
    global imported_image
    if imported_image is not None:
        original_h, original_w = imported_image.shape
        sensor_h, sensor_w = image.shape

        # Compute new dimensions while maintaining aspect ratio
        new_width = int(sensor_w * scale_factor)
        new_height = int(original_h * (new_width / original_w))

        # Ensure the resized image does not exceed sensor dimensions
        new_height = min(new_height, sensor_h)
        new_width = min(new_width, sensor_w)

        resized_image = resize(imported_image, (new_height, new_width), anti_aliasing=True)
        pixel_flux = calculate_image_flux(image_magnitude, sensor_params)
        resized_image *= pixel_flux  # Scale the image based on calculated flux

        # Compute centering positions
        start_x = max(0, (sensor_w - new_width) // 2)
        start_y = max(0, (sensor_h - new_height) // 2)

        # Ensure the final cropped region does not exceed bounds
        end_x = start_x + new_width
        end_y = start_y + new_height

        # Crop to ensure compatibility with NumPy broadcasting
        image[start_y:end_y, start_x:end_x] += resized_image[:end_y - start_y, :end_x - start_x]

    return image

### PSF addition ###

def add_psf(image, x, y, flux, sigma):
    """Add a 2D Gaussian PSF to the image at (x, y) with the given flux."""
    size = int(6 * sigma)
    # Create grid centered on 0
    y_indices, x_indices = np.meshgrid(np.arange(-size//2, size//2+1),
                                       np.arange(-size//2, size//2+1), indexing='ij')
    psf = np.exp(-(x_indices**2 + y_indices**2) / (2 * sigma**2))
    psf /= psf.sum()
    ix, iy = int(y), int(x)
    # Determine sub-image boundaries
    if 0 <= ix < image.shape[0] and 0 <= iy < image.shape[1]:
        x_start, x_end = max(0, ix - size//2), min(image.shape[0], ix + size//2+1)
        y_start, y_end = max(0, iy - size//2), min(image.shape[1], iy + size//2+1)
        sub_psf = psf[:x_end-x_start, :y_end-y_start]
        image[x_start:x_end, y_start:y_end] += flux * sub_psf

### Main image generation function ###

def generate_image(params, binning=False, cosmic_rays=False, sky_background=False, moving_exposures=False, snr_calc=False, pbar=None, win=None, progress_update=None):

    if progress_update: progress_update(10, "Generating empty image...")
    image = np.zeros((int(params["Sensor Height (px)"]), int(params["Sensor Width (px)"])))

    if progress_update: progress_update(20, "Generating random star positions...")
    x_positions, y_positions, magnitudes = [int(params["Sensor Height (px)"])/2], [int(params["Sensor Width (px)"])/2], [params["Star Magnitude"]]

    if progress_update: progress_update(20, "Adding star PSFs to image...")

    if moving_exposures:
        # Divide the exposure into subexposures to simulate camera drift.
        num_steps = 10  # Can be a parameter?

        trail_length = params["Trail Length (px)"]
        drift_angle_rad = np.deg2rad(params["Drift Angle (deg)"])
        dx = trail_length * np.cos(drift_angle_rad) / (num_steps - 1)
        dy = trail_length * np.sin(drift_angle_rad) / (num_steps - 1)
        for step in range(num_steps):
            sub_exposure_time = params["Exposure Time"] / num_steps
            # For each subexposure, add the star signals with a small positional offset
            for x, y, mag in zip(x_positions, y_positions, magnitudes):
                flux = 10 ** (-0.4 * (mag - params["Zero Point"]))
                photons = flux * sub_exposure_time
                electrons = np.random.poisson(photons * params["QE"])
                add_psf(image, x + dx * step, y + dy * step, electrons, params["PSF (sigma)"])
    else:
        # Normal (static) exposure
        for x, y, mag in zip(x_positions, y_positions, magnitudes):
            flux = 10 ** (-0.4 * (mag - params["Zero Point"]))
            photons = flux * params["Exposure Time"]
            electrons = np.random.poisson(photons * params["QE"])
            add_psf(image, x, y, electrons, params["PSF (sigma)"])

    if progress_update:
        progress_update(20, "Adding noise to image...")

    # Add cosmic rays if toggled
    if cosmic_rays:
        image = add_cosmic_rays(image,
                                 num_rays=int(params["Cosmic Ray Count"]),
                                 max_length=int(params["Cosmic Ray Max Length"]),
                                 intensity=int(params["Cosmic Ray Intensity (e-)"]))
    
    # Add sky background if toggled
    if sky_background:
        image = add_sky_background(image,
                                   background_rate=params["Sky Background Rate (e-/px/s)"],
                                   exposure_time=params["Exposure Time"])
    
    # Add dark noise and readout noise
    dark_noise = np.random.poisson(params["Dark Current (e-)"] * params["Exposure Time"], image.shape)
    readout_noise = np.random.normal(params["Readout Noise (e-)"], 1.5, image.shape).astype(int)
    image += dark_noise + readout_noise
    
    if progress_update:
            progress_update(20, "Finalizing image...")

    # Clip the image to sensor's saturation capacity
    image = np.clip(image, 0, params["Saturation Capacity (e-)"]).astype(int)
    
    # Apply 3x3 binning if toggled
    if binning:
        print('Applying 3x3 binning')
        image = apply_binning(image, bin_size=3)
    
    if snr_calc:
        
        print('Calculating SNR')
        # Calculate the signal-to-noise ratio.
        
        if binning: psf = params["PSF (sigma)"] / 3
        else: psf = params["PSF (sigma)"]

        # geometry
        H, W = image.shape
        cx, cy = W//2, H//2
        r_ap = int(3 * psf)           # star aperture
        r_in  = r_ap + 2*max(1, int(psf))   # background annulus
        r_out = r_in + 3*max(1, int(psf))

        yy, xx = np.ogrid[:H, :W]
        rr2 = (xx - cx)**2 + (yy - cy)**2

        ap_mask  = rr2 <= r_ap**2
        bkg_mask = (rr2 >= r_in**2) & (rr2 <= r_out**2)

        # background statistics from annulus
        bkg_vals = image[bkg_mask]
        mu_b  = np.mean(bkg_vals)          # background mean (per pixel)
        sig_b = np.std(bkg_vals, ddof=1)   # background std (per pixel)

        # STAR EXCESS FLUX (sum of background-subtracted pixels in the aperture)
        ap_vals = image[ap_mask]
        N_ap  = ap_vals.size
        N_bkg = bkg_vals.size

        excess_flux = np.sum(ap_vals - mu_b)   # this is the signal (in image units)

        # Empirical noise on that sum: background noise inside aperture + error of background mean
        # Variance of sum from background fluctuations in aperture: N_ap * sig_b^2
        # Variance from background-mean estimate used for subtraction: (N_ap^2) * (sig_b^2 / N_bkg)
        var_emp = N_ap * sig_b**2 + (N_ap**2) * (sig_b**2 / N_bkg)
        noise_emp = np.sqrt(var_emp)

        SNR_empirical = excess_flux / noise_emp

        print("Mask Radius: ", r_ap, ", Signal: ", excess_flux, ", Noise: ", noise_emp, ", SNR: ", SNR_empirical)
        if not np.isnan(SNR_empirical): print(f"SNR: {SNR_empirical:.2f} for Magnitude: {params["Star Magnitude"]:.2f}")

        # This is a simple estimate assuming Poisson noise for CCD/CMOS.
        # signal = params["Exposure Time"] * params["QE"] * 10**(-0.4 * (params["Min Magnitude"] - params["Zero Point"]))
        # noise = np.sqrt(signal + int(sky_background_var.get()) * params["Sky Background Rate (e-/px/s)"] + params["Dark Current (e-)"] * params["Exposure Time"] + params["Readout Noise (e-)"]**2)
        # snr = signal / noise
        # print("Expected Signal: ", signal, "Expected Noise: ", noise, "Expected SNR: ", snr)

    
    return image

### Image saving ###

def save_image(image, filename, format):
    if format == 'png':
        plt.imsave(filename, image, cmap='gray')
    elif format == 'fits':
        fits.writeto(filename, image, overwrite=True)

### GUI Control ###

def run_simulation():
    # Retrieve parameters from the text entries
    win, pbar, status_label = progress_window()

    def update_progress(value, message=None):
        pbar.step(value)
        if message:
            status_label.config(text=message)
            print(message)
        win.update()

    pbar['value'] = 0
    update_progress(0, "Fetching Parameters")
    win.update()

    params = {key: float(entries[key].get()) if entries[key].get() else DEFAULTS[key]
              for key in DEFAULTS}
    params.update({key: float(opt_entries[key].get()) if opt_entries[key].get() else OPTIONAL_PARAMS[key]
              for key in OPTIONAL_PARAMS})

    image = generate_image(params,
                            binning=binning_var.get(),
                            cosmic_rays=cosmic_rays_var.get(),
                            sky_background=sky_background_var.get(),
                            moving_exposures=moving_exposures_var.get(),
                            snr_calc=snr_calc_var.get(),
                            progress_update=update_progress)

    if imgshow_var.get():
        update_progress(9.9, 'Loading image rendering')

        fig = plt.figure(facecolor=(0.2, 0.2, 0.2))
        # plt.style.use('dark_background')
        # fig = plt.figure(figsize=(8, 6))

        ax = fig.add_subplot(111)

        COLOR = 'white'
        mpl.rcParams['text.color'] = COLOR
        mpl.rcParams['axes.labelcolor'] = COLOR
        mpl.rcParams['xtick.color'] = COLOR
        mpl.rcParams['ytick.color'] = COLOR
        mpl.rcParams['axes.edgecolor'] = COLOR
        mpl.rcParams['figure.facecolor'] = (0.2, 0.2, 0.2)
        mpl.rcParams['figure.edgecolor'] = COLOR
        mpl.rcParams['font.family'] = 'monospace'

        ax.set_title("Simulated CMOS Image", color='white')
        ax.set_aspect('equal')

        ax.tick_params(axis='both', left=False, top=False, right=False, bottom=False, labelleft=False, labeltop=False, labelright=False, labelbottom=False)
        spines = ['top', 'bottom', 'left', 'right']
        for spine in spines:
            ax.spines[spine].set_color('white')
            ax.spines[spine].set_linewidth(0.5)

        update_progress(0, "Displaying image...")

        Image = ax.imshow(image, cmap='gray', origin='lower')
        fig.subplots_adjust(bottom=0.2, left=0.1)
        fig.colorbar(Image, ax=ax, label='Electron Count', shrink=0.8)
        vmin = fig.add_axes([0.15, 0.15, 0.7, 0.02])
        vmax = fig.add_axes([0.15, 0.10, 0.7, 0.02])
        vmin_slider = Slider(vmin, 'Min Value', 0, image.max(), valinit=0, color='red')
        vmax_slider = Slider(vmax, 'Max Value', 0, image.max(), valinit=image.max(), color='red')
        
        def update(val):
            new_vmin = vmin_slider.val
            new_vmax = vmax_slider.val
            Image.set_clim(new_vmin, new_vmax)
            fig.canvas.draw_idle() # Redraw the figure
            plt.draw() # Update the canvas

        vmin_slider.on_changed(update)
        vmax_slider.on_changed(update)

        # fig.tight_layout()
        win.destroy()
        plt.show()
    else:
        update_progress(100, "Simulation complete. Image generated.")
        win.destroy()
    # Plot histogram of pixel values
    # mu, sigma = np.mean(image), np.std(image)
    # x = np.linspace(mu - 3*sigma, mu + 3*sigma, 100)
    # plt.plot(x, 6*stats.norm.pdf(x, mu, sigma), linewidth=0.5, label='Gaussian Fit')
    # plt.hist(image.flatten(), bins=100, range=(0, image.max()))
    # plt.yscale('log')
    # plt.show()

def save_file():
    filename = filedialog.asksaveasfilename(defaultextension=".png",
                                            filetypes=[("PNG Image", "*.png"), ("FITS File", "*.fits")])
    if filename:
        format = "fits" if filename.endswith(".fits") else "png"
        params = {key: float(entries[key].get()) if entries[key].get() else DEFAULTS[key]
                  for key in DEFAULTS}
        # Retrieve optional parameters
        params.update({key: float(entries[key].get()) if entries[key].get() else OPTIONAL_PARAMS[key]
                       for key in OPTIONAL_PARAMS})
        image = generate_image(params,
                               binning=binning_var.get(),
                               cosmic_rays=cosmic_rays_var.get(),
                               sky_background=sky_background_var.get(),
                               moving_exposures=moving_exposures_var.get(),
                               snr_calc=snr_calc_var.get())
        save_image(image, filename, format)
        messagebox.showinfo("Success", f"Image saved as {filename}")

def progress_window():
    new_window = tk.Toplevel(root)
    new_window.title("Simulation Progress")

    tk.Label(new_window, text="Simulation in Progress...", font=("TkDefaultFont", 24, 'bold')).grid(row=0, column=0, padx=100, pady=20, sticky='sew')
    status_label = tk.Label(new_window, text="", font=("TkDefaultFont", 20))
    status_label.grid(row=1, column=0, padx=100, sticky='ew')
    # status_label.pack(pady=5)
    
    new_window.geometry("+100+100")
    new_window.resizable(False, False)

    progress_bar = ttk.Progressbar(new_window, mode="determinate", length=300, orient="horizontal", maximum=100)
    progress_bar.grid(row=2, column=0, padx=50, pady=20,sticky='new')
    return new_window, progress_bar, status_label



### GUI Setup ###

root = tk.Tk()
root.title("CMOS Image Simulation GUI")

title_frame = tk.Frame(root)
title_frame.grid(row=0, column=0, columnspan=2, padx=120, pady=10, sticky="nsew")
tk.Label(title_frame, text="CMOS Image Simulation (Single Star)", font=("TkDefaultFont", 24, 'bold')).grid(row=0, column=0, sticky="w")
tk.Label(title_frame, text="by WashU Satellite", font=("TkDefaultFont", 20)).grid(row=1, column=0, sticky="w")
tk.Label(title_frame, text="Version 1.0", font=("TkDefaultFont", 16)).grid(row=2, column=0, sticky="w")

logo_frame = tk.Frame(root)
logo_frame.grid(row=0, column=0, padx=20, pady=10, sticky="nw")
logo_img = Image.open("cmos-sims/NewDLogo.png")
title_frame.update_idletasks()
logo_img = logo_img.resize((title_frame.winfo_height(), title_frame.winfo_height()), Image.LANCZOS)
logo = ImageTk.PhotoImage(logo_img)
logo_label = tk.Label(logo_frame, image=logo)
logo_label.pack()
logo_label.image = logo # reference for garbage collection

left_frame = tk.Frame(root)
left_frame.grid(row=1, column=0, padx=10, sticky="nsew")

# Create parameter entries frame
param_frame = tk.LabelFrame(left_frame, text="Default Parameters", padx=10, pady=10)
param_frame.grid(row=0, column=0, padx=10, sticky="new")

entries = {}
# List default parameters in GUI
for i, (key, value) in enumerate(DEFAULTS.items()):
    tk.Label(param_frame, text=key).grid(row=i, column=0, pady=0, sticky="w")
    entries[key] = tk.Entry(param_frame, width=12)
    entries[key].grid(row=i, column=1)
    entries[key].insert(0, str(value))

button_frame = tk.Frame(root)
button_frame.grid(row=1, column=1, sticky="nsew")

# GUI for optional toggles
toggle_frame = tk.LabelFrame(button_frame, text="Optional Features", padx=10, pady=10)
toggle_frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")

# opt_vals = tk.Frame(toggle_frame).grid(row=0, column=0)
# opt_bools = tk.Frame(toggle_frame).grid(row=0, column=1)

binning_var = tk.BooleanVar(value=False)
cosmic_rays_var = tk.BooleanVar(value=False)
sky_background_var = tk.BooleanVar(value=False)
imgshow_var = tk.BooleanVar(value=True)
moving_exposures_var = tk.BooleanVar(value=False)
snr_calc_var = tk.BooleanVar(value=False)

tk.Checkbutton(toggle_frame, text="Simulate Moving Exposures", variable=moving_exposures_var).grid(row=0, column=0, sticky='w')
tk.Checkbutton(toggle_frame, text="Add Cosmic Rays", variable=cosmic_rays_var).grid(row=1, column=0, sticky='w')
tk.Checkbutton(toggle_frame, text="Add Sky Background", variable=sky_background_var).grid(row=2, column=0, sticky='w')
tk.Checkbutton(toggle_frame, text="Show Image", variable=imgshow_var).grid(row=0, column=1, sticky='w')
tk.Checkbutton(toggle_frame, text="3x3 Binning", variable=binning_var).grid(row=1, column=1, sticky='w')
tk.Checkbutton(toggle_frame, text="Calculate SNR", variable=snr_calc_var).grid(row=2, column=1, sticky='w')

opt_param_frame = tk.LabelFrame(button_frame, text="Optional Parameters", padx=10, pady=10)
opt_param_frame.grid(row=2, column=0, padx=10, sticky="nsew")

opt_entries = {}
for i, (key, value) in enumerate(OPTIONAL_PARAMS.items()):
    tk.Label(opt_param_frame, text=key).grid(row=i, column=0, sticky="w")
    opt_entries[key] = tk.Entry(opt_param_frame, width=12)
    opt_entries[key].grid(row=i, column=1)
    opt_entries[key].insert(0, str(value))

toggle_map = [
    (cosmic_rays_var, ["Cosmic Ray Count", "Cosmic Ray Max Length", "Cosmic Ray Intensity (e-)"]),
    (sky_background_var, ["Sky Background Rate (e-/px/s)"]),
    (moving_exposures_var, ["Trail Length (px)", "Drift Angle (deg)"]),
]

def toggle_entries(var, keys):
    state = "normal" if var.get() else "disabled"
    for key in keys:
        opt_entries[key].config(state=state)

for var, keys in toggle_map:
    var.trace_add("write",
        lambda *_, v=var, ks=keys: toggle_entries(v, ks)
    )
    toggle_entries(var, keys)

for e in opt_entries.values():
    e.config(disabledbackground="#969595",
             disabledforeground="#666666",
             relief="flat",
             bd=1,)

# Create frame for action buttons
button_frame = tk.Frame(root)
button_frame.grid(row=3, column=0, columnspan=2, padx=20, pady=20)
tk.Button(button_frame, text="Run Simulation", command=run_simulation, width=20, height=2).grid(row=0, column=0, padx=5, sticky="nsew")
tk.Button(button_frame, text="Save Image", command=save_file, width=20, height=2).grid(row=0, column=1, padx=5, sticky="nsew")

global current_progress
current_progress = 'Initializing...'

root.mainloop()