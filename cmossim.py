import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from astropy.io import fits

# Constants
PIXEL_SIZE = 3.76e-6*3  # meters
SENSOR_WIDTH = int(9568/3)   # pixels
SENSOR_HEIGHT = int(6380/3)  # pixels
FOCAL_LENGTH = 0.2    # meters
APERTURE = 0.111      # meters
QE = 0.6              # Quantum efficiency (50%)
LAMBDA_NM = 640       # Wavelength in nm
DARK_CURRENT = 0.56   # e-/s
SATURATION_CAPACITY = 51000  # e-
READOUT_NOISE = 1     # e-
FIELD_OF_VIEW = 10    # degrees
MAG_LIMIT = 20        # Dim star magnitude limit
MAG_BRIGHT = 0        # Bright star magnitude limit
PSF_SIGMA = 3         # Gaussian PSF standard deviation in pixels

def mag_limit_calc(D, QE, lambda_nm, pixel_size, sensor_width, sensor_height):
    h = 6.626e-34  # Planck's constant (J·s)
    c = 3.0e8      # Speed of light (m/s)
    lambda_m = lambda_nm * 1e-9  # Convert to meters
    F_Vega = 2.55e-12   # Vega flux in W/m2/Hz
                        # Vega Flux is 2550 Janskys at 640 nm according to:
                        # http://brucegary.net/SED/#:~:text=Since%20a%20zero%20mag%20star%20at%20V,interest%20has%20a%20flux%20of%200.391%20[Jy].&text=The%20last%20column%20is%20used%20to%20calculate,10%2D0.4%20MAG%2C%20where%20MAG%20is%20star%20magnitude.
                        # Vega Flux is 3174 Janskys at 616 nm according to:
                        # https://www.gemini.edu/observing/resources/magnitudes-and-fluxes

    # Aperture area
    A_telescope = np.pi * (D / 2) ** 2

    # Compute photon flux from a mag 0 star
    photon_energy = h * c / lambda_m
    Num_photon_Vega = (F_Vega * 30 * A_telescope * QE) / photon_energy

    # Compute field of view per pixel (square degrees)
    # pixel_fov = (pixel_size / (D * np.pi / 180)) ** 2
    # pixel_fov = (pixel_size / FOCAL_LENGTH) ** 2

    # Compute photon flux per pixel
    # N_photon_Vega_per_pixel = Num_photon_Vega * pixel_fov
    
    # Pixel and sensor areas
    A_pixel = pixel_size ** 2  # Area of one pixel (m^2)
    A_sensor = (sensor_width * pixel_size) * (sensor_height * pixel_size)  # Total sensor area (m^2)

    # Compute photon flux per pixel
    N_photon_Vega_per_pixel = Num_photon_Vega * (A_pixel / A_sensor)

    # Compute corresponding zero-point magnitude for per-pixel scaling
    mag_zero_point_physical = -2.5 * np.log10(N_photon_Vega_per_pixel)
    return mag_zero_point_physical

def simulate_stars(num_stars, mag_bright, mag_limit):
    # Generate random star positions and magnitudes within the field of view
    x_positions = np.random.uniform(0, SENSOR_WIDTH, num_stars)
    y_positions = np.random.uniform(0, SENSOR_HEIGHT, num_stars)
    magnitudes = np.random.uniform(mag_bright, mag_limit, num_stars)
    return x_positions, y_positions, magnitudes

def magnitude_to_flux(magnitude, mag_zero_point=0):
    # Convert magnitude to flux (arbitrary units normalized to the zero point)
    # The zero point is the flux of a mag 0 star in the given aperture and QE
    # Ex: If mag_zero_point = 10, then a magnitude 10 star produces 1 photon per pixel per second.

    flux = 10**(-0.4 * (magnitude - mag_zero_point))
    return flux

def magnitude_to_flux2(magnitude, mag_zero_point=0):
    # Another way to calculate it
    P_vega = 702.0    # photons cm-2 s-1 A-1
                        # https://www.astronomy.ohio-state.edu/martini.10/usefuldata.html
    bandwidth = 6750-6450 # A
    F_vega_area = P_vega * bandwidth # photons cm-2 s-1
    A = np.pi * ((APERTURE * 100)/2)**2 # cm2
    F_vega = F_vega_area * A # photons s-1
    
    flux = F_vega*10**(-0.4 * magnitude)
    return flux

def generate_psf(image, x, y, flux, sigma):
    """
    Add a star with a Gaussian PSF to the image.
    Parameters:
    - image: 2D NumPy array (sensor image)
    - x, y: Pixel coordinates of the star
    - flux: Total photon count to distribute in the PSF
    - sigma: Standard deviation of the Gaussian PSF in pixels
    """
    size = int(6 * sigma)  # Define PSF size (approximately 3σ in each direction)
    # Generate the PSF kernel
    y_indices, x_indices = np.meshgrid(
        np.arange(-size // 2, size // 2 + 1),
        np.arange(-size // 2, size // 2 + 1),
        indexing='ij'
    )
    psf = np.exp(-(x_indices**2 + y_indices**2) / (2 * sigma**2))
    psf /= psf.sum()  # Normalize PSF so total flux matches input
    # Determine star position within the image
    ix, iy = int(y), int(x)
    # Calculate bounds for the image patch where the PSF will be applied
    x_start = max(0, ix - size // 2)
    x_end = min(SENSOR_HEIGHT, ix + size // 2 + 1)
    y_start = max(0, iy - size // 2)
    y_end = min(SENSOR_WIDTH, iy + size // 2 + 1)
    # Ensure the extracted PSF region matches the image region size
    sub_psf_x_start = max(0, size // 2 - ix)
    sub_psf_x_end = sub_psf_x_start + (x_end - x_start)
    sub_psf_y_start = max(0, size // 2 - iy)
    sub_psf_y_end = sub_psf_y_start + (y_end - y_start)
    # Extract the correctly sized sub-region from the PSF
    sub_psf = psf[sub_psf_x_start:sub_psf_x_end, sub_psf_y_start:sub_psf_y_end]
    # Add the Gaussian PSF centered at (ix, iy) while ensuring boundaries are respected
    image[x_start:x_end, y_start:y_end] += flux * sub_psf

def generate_image(exposure_time=1.0, num_stars=1000):
    # Simulate a CMOS image with stars and a Gaussian PSF.

    # Create an empty image
    image = np.zeros((SENSOR_HEIGHT, SENSOR_WIDTH))

    # Simulate star positions and magnitudes
    x_positions, y_positions, magnitudes = simulate_stars(num_stars, MAG_BRIGHT, MAG_LIMIT)

    # Compute the magnitude zero point for the given aperture and QE
    mag_zero_point = mag_limit_calc(APERTURE, QE, LAMBDA_NM, PIXEL_SIZE, SENSOR_WIDTH, SENSOR_HEIGHT)

    # Add stars to the image
    for i, (x, y, mag) in enumerate(zip(x_positions, y_positions, magnitudes)):
        mag_zero_point = 18
        flux = magnitude_to_flux2(mag, mag_zero_point)
        photons = flux * exposure_time
        electrons = np.random.poisson(photons * QE)  # Add photon shot noise
        # electrons = photons * QE
        
        # Print flux and electron values for the first few stars
        if i < 5:  
            print(f"Star {i+1}: Magnitude={mag}, Flux={flux}, Photons={photons}, Electrons={electrons}")

        # Add star with Gaussian PSF
        generate_psf(image, x, y, electrons, PSF_SIGMA)

    # Add dark current and readout noise
    # dark_noise = np.random.normal(DARK_CURRENT * exposure_time, READOUT_NOISE, image.shape)
    dark_noise = np.random.poisson(DARK_CURRENT * exposure_time, image.shape)
    dark_noise += np.random.normal(READOUT_NOISE, 1.5, image.shape).astype(int)
    image += dark_noise

    # Clip at saturation capacity. We aren't simulating blooming yet.
    image = np.clip(image, 0, SATURATION_CAPACITY).astype(int)

    print(type(image[0][0]))
    return image

def save_image(image, filename, format='png'):
    # Save the simulated image to a file
    if format == 'png':
        plt.imsave(filename, image, cmap='gray', format=format)
    elif format == 'fits':
        hdu = fits.PrimaryHDU(image)
        hdu.writeto(filename, overwrite=True)
    else:
        raise ValueError("Unsupported format. Use 'png' or 'fits'.")

# Main simulation
if __name__ == "__main__":
    exposure_time = 10.0  # Default exposure time in seconds.
    num_stars = 1000     # Number of stars to simulate

    # Generate the CMOS image
    cmos_image = generate_image(exposure_time=exposure_time, num_stars=num_stars)

    # Save the image
    # save_image(cmos_image, "cmos_simulation_with_psf.png", format="png")

    print(f"Image max value: {cmos_image.max()}, mean value: {cmos_image.mean()}")
    
    ### Log scale visualization
    # plt.imshow(np.log1p(cmos_image), cmap='gray', origin='lower')  
    # plt.colorbar(label='Log(Electron Count)')
    # plt.title('Simulated CMOS Image with Gaussian PSF')
    # plt.xlabel('Pixel X')
    # plt.ylabel('Pixel Y')
    # plt.show()

    # Linear scale visualization
    plt.imshow(cmos_image, cmap='gray', origin='lower')
    plt.colorbar(label='Electron Count')
    plt.title('Simulated CMOS Image with Gaussian PSF')
    plt.xlabel('Pixel X')
    plt.ylabel('Pixel Y')
    plt.show()