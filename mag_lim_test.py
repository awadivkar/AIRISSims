from cmossim import mag_limit_calc, magnitude_to_flux2
import numpy as np

PIXEL_SIZE = 3.76e-6  # meters
SENSOR_WIDTH = 9568   # pixels
SENSOR_HEIGHT = 6380  # pixels
APERTURE = 0.111      # meters
QE = 0.6              # Quantum efficiency (50%)
LAMBDA_NM = 640       # Wavelength in nm

# mag_lim = mag_limit_calc(APERTURE, QE, LAMBDA_NM, PIXEL_SIZE, SENSOR_WIDTH, SENSOR_HEIGHT)

# print(f'mag_lim_zero: {mag_lim}')
# if mag_lim < 10:
#     print('fuck')
# else:
#     print('we\'re so back')

# def photons_sec_calc(D, QE, lambda_nm, pixel_size, sensor_width, sensor_height, m):
#     P_0 = 995 # Photons/sec*cm^2*A of Vega (not accurate at wavelength)
#     A_telescope = np.pi*(D*100 / 2)**2
#     lambda_width = 6750-6450

print(f'{(magnitude_to_flux2(0))} photons/sec')