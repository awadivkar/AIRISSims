from astroquery.vizier import Vizier
from astropy.coordinates import SkyCoord
import time
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np
from tqdm import tqdm  
import pandas as pd

data = pd.read_csv('cmos-sims/SNR_sampling/mag+expstats.csv')

# plt.figure(figsize=(10, 6))
# plt.scatter(data['Exposure Time'], data['SNR'], alpha=0.1, color='k')
# plt.xlabel('Exposure Time')
# plt.ylabel('SNR')
# plt.show()

plt.hist2d(data['Exposure Time'], data['Magnitude'], bins=100, cmap='plasma')
plt.colorbar(label='Counts')
plt.xlabel('Exposure Time')
plt.ylabel('Magnitude')
plt.title('2D Histogram of Exposure Time vs Magnitude')
plt.show()

print(min(data['Exposure Time']), max(data['Exposure Time']))
print(min(data['Magnitude']), max(data['Magnitude']))

# print(data)