from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import time
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm  

from astropy import units as u


coord = SkyCoord(ra=20*u.deg, dec=20*u.deg, frame='icrs')  # Example coordinates
# radius = 5 * u.deg  # Example radius

fov_deg_h = 6  # Field of view height in degrees
fov_deg_w = 10  # Field of view width in degrees
high_mag_limit = 20  # Example magnitude limit
low_mag_limit = 15  # Example magnitude limit

# vizier = Vizier(columns=['RAJ2000', 'DEJ2000', 'magG'])
# vizier.ROW_LIMIT = -1
# result = vizier.query_region(
#         coord, height=fov_deg_h * u.deg, width=fov_deg_w * u.deg, catalog='igsl3', column_filters={'magG': f'< {high_mag_limit}'})
    
print('Querying Vizier catalog...')
start = time.time()
vizier = Vizier(columns=['RAJ2000', 'DEJ2000', 'magG'])  # Specify the columns you want
vizier.ROW_LIMIT = -1
result = vizier.query_region(coord, height=6 * u.deg, width=10 * u.deg, catalog='igsl3', column_filters={'magG': f'< {high_mag_limit}', 'magG': f'> {low_mag_limit}'})
print('Query completed in {:.2f} seconds'.format(time.time() - start))
print(result[0])
# plt.hist(result[0]['magG'], bins=50, color='blue', alpha=0.7)
for i in tqdm(range(len(result[0]))):
    flux = 10 ** ((result[0]['magG'][i] - 15)/-2.5 + np.log10(5))  # Convert magnitude to flux
    plt.plot(result[0]['RAJ2000'], result[0]['DEJ2000'], 'ro', markersize=result[0]['magG'][i], alpha=0.1)
# plt.yscale('log')
plt.show()
