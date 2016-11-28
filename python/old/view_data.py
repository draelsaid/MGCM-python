from scipy.io import *
from pylab import *

import matplotlib.pyplot as plt
import numpy as np

mgcmdir = "mgcm6"
rundir = "adam_first_run"
mthdir = "m1"
filename = "diagfi.nc"

a = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/%s/%s/%s/%s" % (mgcmdir,rundir,mthdir,filename),'r')

lat = a.variables['lat'][:]
lon = a.variables['lon'][:]
sigma = a.variables['sigma'][:]
time = a.variables['time'][:]
controle = a.variables['controle'][:]
tsurf = a.variables['tsurf'][:]
ps = a.variables['ps'][:]
temp = a.variables['temp'][:]
u = a.variables['u'][:]
v = a.variables['v'][:]
dustq = a.variables['dustq'][:]
dustN = a.variables['dustN'][:]

# Convert sigma to altitude (m)
alt = - 10.8*np.log(sigma)
dustq_s_lat = dustq[:,:,:,0]

# Time series plot of entire mars year, with altitude vs latitude of dust mix ratio to air kg/kg(air) 
fig = plt.figure(1)

m=1
for i in xrange(48):
 ax = plt.subplot(12,4,i+1) 
 im = plt.contourf(lat,alt,dustq_s_lat[i+m,:,:])
 
 ax.yaxis.set_visible(False)
 ax.tick_params(axis='y', which='major', labelsize=6)
 ax.tick_params(axis='x', which='major', labelsize=10)
 
 ax.title('sol %' % (g))
 
 if (i % 4 == 0): # If i is divisible by 4
  ax.yaxis.set_visible(True)
 
 m = m + 2
 if (m==3):
  m=5

fig.subplots_adjust(right=0.8)

cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(im, cax=cbar_ax)
cb.formatter.set_powerlimits((0, 0))
cb.update_ticks()

fig.text(0.5, 0.04, 'Latitude / degrees', ha='center')
fig.text(0.04, 0.5, 'Altitude / m', va='center', rotation='vertical')

fig.suptitle('One month dust to air ratio (every month, sol 0,20,40,60) / kg/kg', fontsize=16)
fig.savefig("plot.jpg", dpi=400)

# Figure of dust IR optical depth (lat vs sol)
fig = plt.figure(2)



