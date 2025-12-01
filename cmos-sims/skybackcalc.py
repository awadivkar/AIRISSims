import numpy as np
import matplotlib.pyplot as plt

def sky_background_electrons_single_rayleigh(
    params,
    solar_elevation_deg,
    bandpass_min=645.0,       # nm
    bandpass_max=675.0,       # nm
    altitude_ft=140000.0,     # feet
    tauR550_sea_level=0.10,   # vertical Rayleigh optical depth at 550 nm (sea level, order-of-mag)
    ozone_tau_above=0.0,      # vertical O3 optical depth above payload in band (Antarctic summer above 42 km: small; set >0 if you want)
    system_throughput=0.85,   # optics*filter transmission (excluding QE)
    view_zenith_deg=0.0       # 0 = zenith view
):
    # --- Constants
    h = 6.62607015e-34
    c = 299792458.0

    # --- Altitude & remaining column
    altitude_m = altitude_ft / 3.28084
    P0 = 1013.25  # hPa
    H  = 8000.0   # m, scale height (approx)
    f_col = np.exp(-altitude_m / H)  # fraction of sea-level column ABOVE the payload

    # --- Band grid (small)
    wl = np.linspace(bandpass_min, bandpass_max, 7) * 1e-9  # m
    wl_nm = wl * 1e9

    # Rayleigh vertical tau above payload at each λ
    tau_s = tauR550_sea_level * f_col * (wl_nm / 550.0)**(-4)

    # Optional ozone (or other absorption) vertical tau above payload
    # Keep it flat across the narrow band unless you have a cross-section curve
    tau_abs = np.full_like(wl, ozone_tau_above)

    # Total vertical tau above payload for attenuation of solar beam
    tau_tot = tau_s + tau_abs

    # --- Solar & viewing geometry
    # Robust solar airmass (Kasten–Young 1989) for the beam; clip at small elevations
    se = np.deg2rad(np.clip(solar_elevation_deg, 0.1, 90.0))
    z_sun = np.pi/2 - se
    mu0 = np.cos(z_sun)
    mu0 = np.clip(mu0, 0.01, 1.0)

    # View zenith (assume zenith by default)
    vz = np.deg2rad(np.clip(view_zenith_deg, 0.0, 75.0))
    mu = np.cos(vz)

    # Scattering angle Θ between Sun and line-of-sight (for zenith viewing, Θ = solar zenith)
    Theta = z_sun if view_zenith_deg == 0.0 else np.arccos(
        np.cos(z_sun)*np.cos(vz) + np.sin(z_sun)*np.sin(vz)*0.0
    )  # relative azimuth=0 assumed; generalize if needed

    # Rayleigh phase function normalization -> 3/(16π) * (1+cos^2Θ)
    phase = (3.0/(16.0*np.pi)) * (1.0 + np.cos(Theta)**2)

    # --- Solar spectrum at TOA (AM0) at ~660 nm
    # Use a reasonable constant across the narrow band; you can swap in a table later.
    # AM0 irradiance around 660 nm is ~1.5–1.9 W/m^2/nm; let’s pick 1.7 W/m^2/nm.
    F0_lambda = 1.7  # W/m^2/nm
    F0 = np.full_like(wl_nm, F0_lambda)  # W/m^2/nm

    # Attenuate the incoming solar beam by the small column above payload
    trans_solar = np.exp(-tau_tot / mu0)  # elementwise

    # Single-scattered spectral radiance (W m^-2 nm^-1 sr^-1)
    # For small tau_s, (1 - e^{-tau_s/mu})/mu ≈ tau_s/mu
    L_lambda = F0 * trans_solar * phase * (tau_s / mu)

    # Convert to photon radiance and integrate over band
    E_ph = (h * c) / wl   # J per photon
    photon_radiance_per_nm = L_lambda / E_ph  # photons m^-2 s^-1 sr^-1 nm^-1

    # Integrate over the passband (nm)
    dlam_nm = (bandpass_max - bandpass_min) / (len(wl_nm)-1)
    photon_radiance = np.trapz(photon_radiance_per_nm, wl_nm)  # photons m^-2 s^-1 sr^-1

    # --- Pixel solid angle
    px = params["Pixel Size (um)"] * 1e-6       # m
    sensor_w = params["Sensor Width (px)"] * px # m
    fov_rad = np.deg2rad(params["Field of View (deg)"])  # full horizontal FOV
    f = sensor_w / fov_rad                      # effective focal length (small-angle)
    pixel_solid_angle = (px / f)**2             # sr

    # --- Effective collecting area and throughput
    D = params["Aperture"]                      # m
    A = np.pi * (D/2.0)**2                      # m^2
    T = system_throughput                       # optics*filter (QE handled below)

    # --- Electrons per pixel
    texp = params["Exposure Time"]              # s
    QE   = params["QE"]                         # detector QE (0..1)

    photons_per_pixel = photon_radiance * A * pixel_solid_angle * texp * T
    electrons_per_pixel = photons_per_pixel * QE

    return float(electrons_per_pixel)

# --- Example Usage ---

params = {
    "Wavelength (nm)": 660,           # Not directly used here, but could be informative
    "Exposure Time": 1,              # seconds
    "QE": 0.6,                        # Quantum efficiency
    "Aperture": 0.111,                # Aperture diameter in meters (example value)
    "Pixel Size (um)": 3.76,           # in microns
    "Sensor Width (px)": 9568,         # number of pixels across sensor width
    "Field of View (deg)": 10         # degrees
}

# For a solar elevation of 10° (the Sun is low on the horizon)
# electrons = calculate_sky_background_electrons(params, solar_elevation_deg=10)
# print("Estimated sky background (electrons per pixel):", electrons)

x = np.linspace(0, 90, 100)
y = [sky_background_electrons_single_rayleigh(params, solar_elevation_deg=el) for el in x]
plt.plot(x, y)
# plt.yscale('log')
plt.xlabel('Solar Elevation (degrees)')
plt.ylabel('Sky Background (electrons per pixel)')
plt.title('Sky Background vs Solar Elevation at 140,000 ft')
plt.grid(True)
plt.show()
