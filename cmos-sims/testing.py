from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import time
import matplotlib.pyplot as plt
import numpy as np
from tqdm import tqdm  

from astropy import units as u

vals = []
for x in tqdm(range(1000)): vals.append(np.random.poisson(100))

plt.hist(vals)
plt.show()