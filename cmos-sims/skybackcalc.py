import numpy as np
import matplotlib.pyplot as plt

def calculate_sky_background_electrons(params,
                                         solar_elevation_deg,
                                         bandpass_min=645,    # nm
                                         bandpass_max=675,    # nm
                                         altitude_ft=140000,  # Camera altitude in feet
                                         tau0=0.1,            # Reference optical depth at 550 nm
                                         ref_wavelength=550   # nm, reference wavelength
                                        ):
    # 1. Convert the altitude from feet to meters.
    altitude_m = altitude_ft / 3.281  # ~42,700 m
    
    # 2. Estimate atmospheric pressure at altitude using an exponential model.
    P0 = 1013  # hPa at sea level
    scale_height = 8000  # meters (typical scale height)
    P_alt = P0 * np.exp(-altitude_m / scale_height)
    
    # Vertical column factor (fraction of atmospheric column remaining above the camera)
    vertical_column_factor = P_alt / P0  # e.g., ~0.005 for 42,700 m
    
    # 3. Effective airmass for the incoming solar light.
    solar_elev_rad = np.deg2rad(solar_elevation_deg)
    effective_airmass = vertical_column_factor / np.sin(solar_elev_rad)
    
    # 4. Bandpass details.
    delta_lambda = bandpass_max - bandpass_min  # nm, should be 30 nm
    central_wavelength = (bandpass_max + bandpass_min) / 2.0  # nm, ~660 nm
    
    # 5. Solar spectral irradiance at central wavelength (W/m²/nm)
    # You may update this value based on a more precise solar spectrum.
    F_lambda = 1.9  # W/m²/nm (example value at 660 nm)
    
    # 6. Wavelength-dependent Rayleigh optical depth.
    tau_rayleigh = tau0 * effective_airmass * (central_wavelength / ref_wavelength)**(-4)
    
    # 7. Scattered radiance into the line of sight (W/m²/nm/sr).
    # (The factor 1/(4*pi) comes from assuming roughly isotropic scattering.)
    radiance = (F_lambda * tau_rayleigh) / (4 * np.pi)
    
    # 8. Convert radiance to photon radiance.
    h = 6.626e-34  # Planck's constant (J s)
    c = 3e8        # Speed of light (m/s)
    central_wavelength_m = central_wavelength * 1e-9  # convert nm to m
    energy_per_photon = h * c / central_wavelength_m  # Joules per photon
    photon_radiance = radiance / energy_per_photon  # photons/m²/s/nm/sr
    
    # Integrate over the bandpass.
    photon_radiance_total = photon_radiance * delta_lambda  # photons/m²/s/sr
    
    # 9. Determine the pixel's solid angle.
    pixel_size = params["Pixel Size (um)"] * 1e-6  # Convert from µm to m
    sensor_width_px = params["Sensor Width (px)"]
    sensor_width_m = sensor_width_px * pixel_size
    fov_rad = np.deg2rad(params["Field of View (deg)"])
    
    # Estimate the effective focal length from the sensor width and field of view.
    focal_length = sensor_width_m / fov_rad  
    # Angular size of one pixel (in radians).
    pixel_angle = pixel_size / focal_length  
    # Solid angle per pixel (in steradians).
    pixel_solid_angle = pixel_angle**2
    
    # 10. Telescope effective area.
    # The "Aperture" parameter is taken to be the diameter in meters.
    aperture_diameter = params["Aperture"]
    area = np.pi * (aperture_diameter / 2)**2  # in m²
    
    # 11. Exposure time and sensor QE.
    exposure_time = params["Exposure Time"]  # seconds
    QE = params["QE"]
    
    # 12. Calculate the number of photons per pixel.
    photons_per_pixel = photon_radiance_total * area * pixel_solid_angle * exposure_time
    
    # 13. Convert to electrons using the sensor's quantum efficiency.
    electrons_per_pixel = photons_per_pixel * QE
    
    return electrons_per_pixel

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
electrons = calculate_sky_background_electrons(params, solar_elevation_deg=10)
print("Estimated sky background (electrons per pixel):", electrons)
# x = np.linspace(0, 90, 100)
# y = [calculate_sky_background_electrons(params, solar_elevation_deg=el) for el in x]
# plt.plot(x, y)
# plt.xlabel("Solar Elevation (degrees)")
# plt.ylabel("Sky Background (electrons per pixel)")
# plt.title("Sky Background vs. Solar Elevation")
# plt.grid(True)
# plt.show()   