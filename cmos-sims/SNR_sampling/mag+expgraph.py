import pandas as pd
import numpy as np
from matplotlib.colors import LogNorm, Normalize, SymLogNorm
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def snr_model(X, a, b, c):
    t_exp, mag = X
    return a * (t_exp ** b) * (10 ** (-0.4 * c * mag))

# Load the CSV file
file_path = 'cmos-sims/SNR_sampling/mag+expstats.csv'
data = pd.read_csv(file_path)

# Extract columns
exposure = data['Exposure Time']
magnitude = data['Magnitude']
snr = data['SNR']

mask = np.isfinite(exposure) & np.isfinite(magnitude) & np.isfinite(snr) & (snr > 0)
popt, pcov = curve_fit(snr_model, (exposure[mask], magnitude[mask]), snr[mask],
                       p0=[8000, 0.8, 1])
a_fit, b_fit, c_fit = popt
print(f"Best-fit parameters:\n  a = {a_fit:.5e}\n  b = {b_fit:.5e}\n  c = {c_fit:.5e}")

# Define the bins for Exposure and Magnitude
n_bins = 80
exposure_bins = np.linspace(exposure.min(), exposure.max(), n_bins)  # Adjust the number of bins as needed
magnitude_bins = np.linspace(magnitude.min(), magnitude.max(), n_bins)

# Create a 2D histogram for SNR
heatmap, xedges, yedges = np.histogram2d(exposure, magnitude, bins=[exposure_bins, magnitude_bins], weights=snr)

# Normalize by number of samples per bin
counts, _, _ = np.histogram2d(exposure, magnitude, bins=[exposure_bins, magnitude_bins])
heatmap = np.divide(heatmap, counts, out=np.zeros_like(heatmap), where=counts != 0)

# Mask zeros to avoid log problems
# heatmap[heatmap <= 0] = heatmap[heatmap > 0].min()
heatmap[heatmap == 0] = np.nan  # Mask zeros for better visualization

fig, ax = plt.subplots(2, 4, figsize=(20, 10))

img = ax[0][0].imshow(
    heatmap.T,
    origin='lower',
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='viridis',
    norm=Normalize(vmin=1, vmax=np.nanmax(heatmap))
)

# --- Add contour at SNR = 5 ---
X, Y = np.meshgrid(exposure_bins[:-1], magnitude_bins[:-1])
levels = [3, 5, 10]
contour = ax[0][0].contour(
    X, Y, heatmap.T, levels=levels,
    colors='red', linewidths=1.5, linestyles='--'
)

ax[0][0].clabel(contour, fmt={3: 'SNR=3', 5: 'SNR=5', 10: 'SNR=10'}, colors='red', fontsize=10)

fig.colorbar(img, norm=Normalize(vmin=1, vmax=np.nanmax(heatmap)), label='SNR (log scale)', ax=ax[0][0])
ax[0][0].set_xlabel('Exposure Time')
ax[0][0].set_ylabel('Magnitude')
ax[0][0].set_title('SNR Heatmap')

img_log = ax[1][0].imshow(
    heatmap.T,
    origin='lower',
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='viridis',
    norm=LogNorm(vmin=heatmap[heatmap > 0].min(), vmax=np.nanmax(heatmap))
)
# --- Add contour at SNR = 5 ---
contour_log = ax[1][0].contour(
    X, Y, heatmap.T, levels=levels,
    colors='red', linewidths=1.5, linestyles='--'
)
ax[1][0].clabel(contour_log, fmt={3: 'SNR=3', 5: 'SNR=5', 10: 'SNR=10'}, colors='red', fontsize=10)
fig.colorbar(img_log, label='SNR (log scale)', ax=ax[1][0])
ax[1][0].set_xlabel('Exposure Time')
ax[1][0].set_ylabel('Magnitude')
ax[1][0].set_title('SNR Heatmap (Log Scale)')

# 1) Build centers (length n)
xcenters = 0.5*(exposure_bins[:-1] + exposure_bins[1:])
ycenters = 0.5*(magnitude_bins[:-1] + magnitude_bins[1:])

# 2) Grid on centers
Xc, Yc = np.meshgrid(xcenters, ycenters, indexing='xy')

# 3) Evaluate model on the same center grid
Zg = snr_model((Xc, Yc), *popt)        # Zg.shape == Xc.shape == Yc.shape
print(Xc.shape, Yc.shape, Zg.shape)

model = ax[0][1].imshow(
    Zg,
    origin='lower',
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='viridis',
    norm=Normalize(vmin=Zg[Zg > 0].min(), vmax=Zg.max())
)
model_cont = ax[0][1].contour(Xc, Yc, Zg, levels=levels, colors='red',
              linewidths=1.5, linestyles='--')
ax[0][1].clabel(model_cont, fmt={3: 'SNR=3', 5: 'SNR=5', 10: 'SNR=10'}, colors='red', fontsize=10)
fig.colorbar(model, label='Predicted SNR', ax=ax[0][1])
ax[0][1].set_xlabel('Exposure Time')
ax[0][1].set_ylabel('Magnitude')
ax[0][1].set_title('Fitted SNR Model and SNR=5 Threshold')

model_log = ax[1][1].imshow(
    Zg,
    origin='lower',
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='viridis',
    norm=LogNorm(vmin=Zg[Zg > 0].min(), vmax=Zg.max())
)
model_log_cont = ax[1][1].contour(Xc, Yc, Zg, levels=levels, colors='red',
              linewidths=1.5, linestyles='--')
ax[1][1].clabel(model_log_cont, fmt={3: 'SNR=3', 5: 'SNR=5', 10: 'SNR=10'}, colors='red', fontsize=10)
fig.colorbar(model_log, label='Predicted SNR', ax=ax[1][1])
ax[1][1].set_xlabel('Exposure Time')
ax[1][1].set_ylabel('Magnitude')
ax[1][1].set_title('Fitted SNR Model and SNR=5 Threshold (Log Scale)')

residuals = (heatmap - Zg.T)
residuals[~np.isfinite(residuals)] = 0  # Mask NaNs
residuals[heatmap == 0] = 0  # Mask bins with no data
res_img = ax[0][2].imshow(
    residuals.T,
    origin='lower',
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='bwr',
    vmin=-np.nanmax(np.abs(residuals)),
    vmax=np.nanmax(np.abs(residuals))
)
fig.colorbar(res_img, label='Residuals', ax=ax[0][2])
ax[0][2].set_xlabel('Exposure Time')
ax[0][2].set_ylabel('Magnitude')
ax[0][2].set_title('Residuals (Data - Model)') 

res_img_log = ax[1][2].imshow(
    residuals.T,
    origin='lower', 
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='bwr',
    norm=SymLogNorm(linthresh=1, linscale=1, vmin=-np.nanmax(np.abs(residuals)), vmax=np.nanmax(np.abs(residuals)))
)
fig.colorbar(res_img_log, label='Residuals', ax=ax[1][2])
ax[1][2].set_xlabel('Exposure Time')
ax[1][2].set_ylabel('Magnitude')
ax[1][2].set_title('Residuals (Data - Model) (Log Scale)')

relative_residuals = residuals / Zg.T
relative_residuals[~np.isfinite(relative_residuals)] = 0  # Mask NaNs
relative_residuals[heatmap == 0] = 0  # Mask bins with no data
res_rel_img = ax[0][3].imshow(
    relative_residuals.T,
    origin='lower',
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='bwr',
    vmin=-np.nanmax(np.abs(relative_residuals)),
    vmax=np.nanmax(np.abs(relative_residuals))
)
fig.colorbar(res_rel_img, label='Residuals', ax=ax[0][3])
ax[0][3].set_xlabel('Exposure Time')
ax[0][3].set_ylabel('Magnitude')
ax[0][3].set_title('Relative Residuals (Data - Model)')

res_rel_img_log = ax[1][3].imshow(
    relative_residuals.T,
    origin='lower', 
    aspect='auto',
    extent=[exposure_bins[0], exposure_bins[-1], magnitude_bins[0], magnitude_bins[-1]],
    cmap='bwr',
    norm=SymLogNorm(linthresh=0.1, linscale=1, vmin=-np.nanmax(np.abs(relative_residuals)), vmax=np.nanmax(np.abs(relative_residuals)))
)
fig.colorbar(res_rel_img_log, label='Relative Residuals', ax=ax[1][3])
ax[1][3].set_xlabel('Exposure Time')
ax[1][3].set_ylabel('Magnitude')
ax[1][3].set_title('Relative Residuals (Data - Model) (Log Scale)')

plt.tight_layout()
plt.savefig('cmos-sims/SNR_sampling/mag_exp_heatmap_with_fits.png', dpi=300)
