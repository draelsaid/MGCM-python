# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016

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

# Use tex font
#rc('text',usetex=True)
# Change all fonts to 'Computer Modern'
#rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})

# Initialise dictionaries - due to data size
Ls_m = {}
psa, psb = {}, {}
presa, presb = {}, {} 
tempa, tempb = {}, {} 
vb, ub = {}, {}
dustqa, dustqb = {}, {}
dustNa, dustNb = {}, {}
rhoa, rhob = {}, {}
fluxsurflwa, fluxsurflwb = {}, {}
fluxsurfswa, fluxsurfswb = {}, {}
fluxtoplwa, fluxtoplwb = {}, {}
fluxtopswa, fluxtopswb = {}, {}
taua, taub = {}, {} 
rdusta, rdustb = {}, {} 
lw_htrta, lw_htrtb = {}, {} 
sw_htrta, sw_htrtb = {}, {} 

# Number of months in comparison (always add 1 because of Python indexing)
Months = 13

# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,Months):
 mgcm = "MGCM_v5-1"
 rundirb = "a_ref4"
 month = "m%i" % (i)
 filename = "diagfi.nc"

 b = netcdf.netcdf_file("/padata/alpha/users/aes442/RUNS/R-%s/%s/%s/%s" % (mgcm,rundirb,month,filename),'r')

 lat = b.variables['lat'][:]
 lon = b.variables['lon'][:]
 sigma = b.variables['sigma'][:]
 t_m = b.variables['time'][:]
 Ls_m[i] = b.variables['Ls'][:]
 
 psb[i] = b.variables['ps'][:]
 presb[i]= b.variables['pressure'][:]
 tempb[i] = b.variables['temp'][:]
 ub[i] = b.variables['u'][:]
 vb[i] = b.variables['v'][:]
 dustqb[i] = b.variables['dustq'][:]
 dustNb[i] = b.variables['dustN'][:]
 rhob[i] = b.variables['rho'][:]
 fluxsurflwb[i] = b.variables['fluxsurf_lw'][:]
 fluxsurfswb[i] = b.variables['fluxsurf_sw'][:]
 fluxtoplwb[i] = b.variables['fluxtop_lw'][:]
 fluxtopswb[i] = b.variables['fluxtop_sw'][:]
 taub[i] = b.variables['taudustvis'][:]
 rdustb[i] = b.variables['reffdust'][:]
 lw_htrtb[i] = b.variables['lw_htrt'][:]
 sw_htrtb[i] = b.variables['sw_htrt'][:]
 
print ("Latitude: %i ||" % (lat.shape)), ("Longitude: %i ||" % (lon.shape)), ("Model levels: %i ||" % (sigma.shape))

# Get time dimension length
n = 0
for i in xrange(1,len(psb)+1,1): # len(psa) gives the number of months
 n = n + len(dustqb[i])          # len(dustqa[i]) gives the number of time steps in each month. Different variable used as a cross-check of dimension consistency.
print ("Total time steps: %i" % (n))
### Create new time dimensions

## Sols vector
mth_s = [61, 66, 66, 65, 60, 54, 50, 46, 47, 47, 51, 46] # Sols per Mars month
sol_s = 0
for i in xrange(0,Months-1):
 sol_s = sol_s + mth_s[i]
t = np.linspace(1,sol_s,n)

## Ls vector
Ls_s = (Months-1)*30 # Number of solar longitudes
Ls = np.zeros((n))

# Method 2 grabs Ls's from model (has bugs, but can be ironed out)
p=0
for i in xrange(1,len(Ls_m)+1,1):
 gg = Ls_m[i]
 for j in xrange(gg.shape[0]):
  Ls[p] = gg[j]
  p = p + 1
Ls = np.roll(Ls,5)
Ls[-1] = np.ceil(Ls[-2])
Ls[:6] = np.linspace(np.floor(Ls[5]),Ls[5],6)
print Ls[:8], Ls[-8:]

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

#data
lvl = 0
latt1 = 17
latt2 = 0
lonn = 36

y = alt

# Common axis labels
cmap = mpl.cm.Accent

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(ub[j].shape[0]):
  ax = axr.plot(ub[j][k,:,latt1,lonn], y, alpha=0.15, linewidth=1.5, color=cmap(1))
  
plt.axis([-400, 250, 0, 100])
axr.set_xlabel('Zonal wind velocity / m/s', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('u_profile.png')

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(vb[j].shape[0]):
  ax = axr.plot(vb[j][k,:,latt1,lonn], y, alpha=0.15, linewidth=1.5, color=cmap(1))
  
plt.axis([-250, 250, 0, 100])
axr.set_xlabel('Meridional wind velocity / m/s', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('v_profile.png')

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(presb[j].shape[0]):
  ax = axr.plot(presb[j][k,:,latt1,lonn], y, alpha=0.15, linewidth=1.5, color=cmap(1))
  
plt.axis([0, 750, 0, 100])
axr.set_xlabel('Pressure / Pa', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('pressure_profile.png')

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(rhob[j].shape[0]):
  ax = axr.plot(rhob[j][k,:,latt1,lonn], y, alpha=0.15, linewidth=1.5, color=cmap(1))
  
plt.axis([0, 0.02, 0, 100])
axr.set_xlabel('Density / kg/$m^3$', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('density_profile.png')

