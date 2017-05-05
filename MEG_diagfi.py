# Pulls out a year's worth of data from the MGCM diagfi.nc for MEG validation
# Adam El-Said OU 04/2017

import matplotlib as mpl
#mpl.use('Agg') # removes need for X-Server (sshing graphics in linux). For qsub only.

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

from mars_time import MarsTime
from scipy.io import *
from matplotlib import cm,ticker
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import *
from scipy.io import *

# Initialise dictionaries - due to data size
Ls_m = {}
psb = {}
presbc, presb = {}, {} 
tempb = {}
vb, ub = {}, {}
rhob = {}
taub = {}

Months = 13

for i in xrange(1,Months):
 mgcm = "MGCM_v5-1"
 rundirb = "T31_9_Marsyears/MY28"
 month = "m%i" % (i)
 filename = "diagfi.nc"
 
 b = netcdf.netcdf_file("/padata/alpha/users/aes442/RUNS/R-%s/%s/%s/%s" % (mgcm,rundirb,month,filename),'r')

 lat = b.variables['lat'][:]
 lon = b.variables['lon'][:]
 sigma = b.variables['sigma'][:]
 t_m = b.variables['time'][:]
# Ls_m[i] = b.variables['Ls'][:]

 psb[i] = b.variables['ps'][:]
 presb[i]= b.variables['pressure'][:]
 tempb[i] = b.variables['temp'][:]
 ub[i] = b.variables['u'][:]
 vb[i] = b.variables['v'][:]
 rhob[i] = b.variables['rho'][:]

# for a run without pressure
 presbc[i] = tempb[i]*0.
 for k in xrange(sigma.shape[0]):
  presbc[i][:,k,:,:] = psb[i][:,:,:]*sigma[k]
 
print ("Latitude: %i || Longitude: %i || Model levels: %i || Months: %i || " % (lat.shape[0], lon.shape[0], sigma.shape[0], Months-1))

# Get time dimension length
n = 0
for i in xrange(1,len(psb)+1,1): # len(psa) gives the number of months
 n = n + len(tempb[i])          # len(dustqa[i]) gives the number of time steps in each month. Different variable used as a cross-check of dimension consistency.
print ("Total time steps: %i" % (n))
### Create new time dimensions

## Sols vector
mth_s = [61, 66, 66, 65, 60, 54, 50, 46, 47, 47, 51, 46] # Sols per Mars month
sol_s = 0
for i in xrange(0,Months-1):
 sol_s = sol_s + mth_s[i]
t = np.linspace(1,sol_s,n)

# Calculate approximate HEIGHT from sigma (km)
alt = np.zeros((sigma.shape[0]))
for i in xrange(len(sigma)):
 alt[i] = -10.8*np.log(sigma[i])

print "PLOTTING....."

# lat = 87.49999, 82.49999, 77.5, 72.5, 67.5, 62.5, 57.5, 52.5, 47.5, 42.5,
#    37.5, 32.5, 27.5, 22.5, 17.5, 12.5, 7.500001, 2.500001, -2.500001,
#    -7.500003, -12.5, -17.5, -22.5, -27.5, -32.5, -37.5, -42.5, -47.5, -52.5,
#    -57.5, -62.5, -67.5, -72.5, -77.5, -82.49999, -87.49999 ;

# lon = -180, -175, -170, -165, -160, -155, -150, -145, -140, -135, -130,
#    -125, -120, -115, -110, -105, -100, -95, -90, -84.99999, -80, -75, -70,
#    -65, -60, -55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10,
#    15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 84.99999, 90, 95,
#    100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165,
#    170, 175 ;

## Data
lt = 45.
ln = 90.

latt = lat[lat>lt].shape[0]
lonn = -lon[lon>ln].shape[0]-1

y = alt

# Common axis labels
cmap = mpl.cm.Accent

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(ub[j].shape[0]):
  ax = axr.plot(ub[j][k,:,latt,lonn], y, alpha=0.05, linewidth=1.5, color=cmap(1))
  
plt.axis([np.min(ub[j][:,:,latt,lonn])-10, np.ceil(np.max(ub[j][:,:,latt,lonn])/10.)*10., 0, 100])
axr.set_xlabel('Zonal wind velocity / m/s', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('MGCM_%s_%i-%i_uprof.png' % (rundirb[-4:], lat[latt], lon[lonn]))

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(vb[j].shape[0]):
  ax = axr.plot(vb[j][k,:,latt,lonn], y, alpha=0.05, linewidth=1.5, color=cmap(1))
  
plt.axis([np.min(vb[j][:,:,latt,lonn])-10, np.ceil(np.max(vb[j][:,:,latt,lonn])/10.)*10., 0, 100])
axr.set_xlabel('Meridional wind velocity / m/s', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('MGCM_%s_%i-%i_vprof.png' % (rundirb[-4:], lat[latt], lon[lonn]))

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(presb[j].shape[0]):
  ax = axr.plot(presb[j][k,:,latt,lonn], y, alpha=0.05, linewidth=1.5, color=cmap(1))
  
plt.axis([0, np.ceil(np.max(presb[j][:,:,latt,lonn])/10.)*10., 0, 100])
axr.set_xlabel('Pressure / Pa', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('MGCM_%s_%i-%i_presprof.png' % (rundirb[-4:], lat[latt], lon[lonn]))

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(rhob[j].shape[0]):
  ax = axr.plot(rhob[j][k,:,latt,lonn], y, alpha=0.05, linewidth=1.5, color=cmap(1))
  
plt.axis([0, np.max(rhob[j][:,:,latt,lonn])+0.005, 0, 100])
axr.set_xlabel('Density / kg/$m^3$', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('MGCM_%s_%i-%i_densprof.png' % (rundirb[-4:], lat[latt], lon[lonn]))

