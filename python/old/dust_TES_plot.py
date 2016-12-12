import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import math

from mars_time import MarsTime
from scipy.io import *
from matplotlib import cm,ticker
from plt_timeseries import plt_timeseries, MidpointNormalize

# Moving average 
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n
    
# Import data from Luca's TES dust files for comparison
a = netcdf.netcdf_file('/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/datafile/dust_MY24.nc','r')

lat_s = a.variables['latitude'][:]
lon_s = a.variables['longitude'][:]
t   = a.variables['Time'][:]
d   = a.variables['dustop'][:]

lat = np.linspace(-90,90,lat_s.shape[0])
lon = np.linspace(-180,180,lon_s.shape[0])

# Zonal averaging (collapses d to 2 dimensions)
d_z = d.sum(axis=2)/d.shape[2] 
 
n=10
d_avg=np.zeros((d_z.shape[0]-(n-1),d_z.shape[1]))
for i in xrange(0,lat_s.shape[0]):
 d_avg[:,i] = moving_average(d_z[:,i],n)

Ls_avg = np.linspace(0,360,t.shape[0]-(n-1))

## PLOT
data = np.matrix.transpose(d_avg)
fig = plt.figure()
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(Ls_avg,lat,data, norm=MidpointNormalize(midpoint=0.), cmap='PuRd')
plt.xlabel('Solar longitude / degrees',fontsize=18)
plt.ylabel('Latitude / degrees',fontsize=18)

# Colour bar
cb = plt.colorbar(format='%.2f')
cb.set_label('Dust optical depth / SI')
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator

plt.axis('tight')
plt.savefig('dust_MY24.png')
