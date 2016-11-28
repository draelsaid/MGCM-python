# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016

from scipy.io import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm,ticker
from matplotlib.ticker import LinearLocator

import matplotlib.pyplot as plt
import numpy as np
import time

from shiftedColorMap import *

# Moving average 
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Initialise dictionaries - due to data size
psa, psb = {}, {} 
tempa, tempb = {}, {}
ua, ub = {}, {}
va, vb = {}, {}
dustqa, dustqb = {}, {}
dustNa, dustNb = {}, {}
rhoa, rhob = {}, {}
fluxsurflwa, fluxsurflwb = {}, {}
fluxsurfswa, fluxsurfswb = {}, {}
fluxtoplwa, fluxtoplwb = {}, {}
fluxtopswa, fluxtopswb = {}, {}
 
# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,13,1):
 rundira = "a_ds"
 rundirb = "a_benchmark" 
 month = "m%i" % (i)
 filename = "diagfi.nc"
 
 a = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/runs/%s/%s/%s" % (rundira,month,filename),'r')
 b = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/runs/%s/%s/%s" % (rundirb,month,filename),'r')

 lat = a.variables['lat'][:]
 lon = a.variables['lon'][:]
 sigma = a.variables['sigma'][:]
 t_m = a.variables['time'][:]
# Ls = a.variables['Ls'][:]
 
 psa[i] = a.variables['ps'][:]
 tempa[i] = a.variables['temp'][:]
 ua[i] = a.variables['u'][:]
 va[i] = a.variables['v'][:]
 dustqa[i] = a.variables['dustq'][:]
 dustNa[i] = a.variables['dustN'][:]
 rhoa[i] = a.variables['rho'][:]
 fluxsurflwa[i] = a.variables['fluxsurf_lw'][:]
 fluxsurfswa[i] = a.variables['fluxsurf_sw'][:]
 fluxtoplwa[i] = a.variables['fluxtop_lw'][:]
 fluxtopswa[i] = a.variables['fluxtop_sw'][:]

 psb[i] = b.variables['ps'][:]
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

# Get time dimension length
n = 0
for i in xrange(1,len(psa)+1,1): # len(psa) gives the number of months
 n = n + len(dustqa[i])          # len(dustqa[i]) gives the number of time steps in each month. Different variable used as a cross-check of dimension consistency.

print n,t_m.shape[0],len(psa),len(dustqa[1])

t = np.arange(1,669,0.25) # Create new time dimension variable. Be careful to change interval if writediagfi interval is changed in runscript

# Create all other variables, with altitude dimension removed, to save time
ps_a, ps_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
temp_a, temp_b   = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
u_a, u_b         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
v_a, v_b         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustq_a, dustq_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustN_a, dustN_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
rho_a, rho_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_lwa, fluxsurf_lwb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_swa, fluxsurf_swb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_lwa, fluxtop_lwb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_swa, fluxtop_swb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))

m=0
for k in xrange(1,len(psa)+1,1):
# Dummy variables - extract month k from dictionary
 dm1, dm1b = psa[k], psb[k]
 dm2, dm2b = tempa[k], tempb[k]
 dm3, dm3b = ua[k], ub[k]
 dm4, dm4b = va[k], vb[k]
 dm5, dm5b = dustqa[k], dustqb[k]
 dm6, dm6b = dustNa[k], dustNb[k]
 dm7, dm7b = rhoa[k], rhob[k]
 dm8, dm8b = fluxsurflwa[k], fluxsurflwb[k]
 dm9, dm9b = fluxsurfswa[k], fluxsurfswb[k]
 dm10, dm10b = fluxtoplwa[k], fluxtoplwb[k]
 dm11, dm11b = fluxtopswa[k], fluxtopswb[k]
 
 # Dummy variables - reduce 4D data to 3D data by selecting surface only (time,lat,lon) 
 dmm2, dmm2b = dm2[:,1,:,:], dm2b[:,1,:,:]
 dmm3, dmm3b = dm3[:,1,:,:], dm3b[:,1,:,:]
 dmm4, dmm4b = dm4[:,1,:,:], dm4b[:,1,:,:]
 dmm5, dmm5b = dm5[:,1,:,:], dm5b[:,1,:,:]
 dmm6, dmm6b = dm6[:,1,:,:], dm6b[:,1,:,:]
 dmm7, dmm7b = dm7[:,1,:,:], dm7b[:,1,:,:] 
 
 for i in xrange(dm2.shape[0]):
  ps_a[m,:,:], ps_b[m,:,:]                = dm1[i,:,:], dm1b[i,:,:]
  temp_a[m,:,:], temp_b[m,:,:]            = dmm2[i,:,:], dmm2b[i,:,:]
  u_a[m,:,:], u_b[m,:,:]                  = dmm3[i,:,:], dmm3b[i,:,:]
  v_a[m,:,:], v_b[m,:,:]                  = dmm4[i,:,:], dmm4b[i,:,:]
  dustq_a[m,:,:], dustq_b[m,:,:]          = dmm5[i,:,:], dmm5b[i,:,:]
  dustN_a[m,:,:], dustN_b[m,:,:]          = dmm6[i,:,:], dmm6b[i,:,:]
  rho_a[m,:,:], rho_b[m,:,:]              = dmm7[i,:,:], dmm7b[i,:,:]
  fluxsurf_lwa[m,:,:], fluxsurf_lwb[m,:,:]= dm8[i,:,:], dm8b[i,:,:]
  fluxsurf_swa[m,:,:], fluxsurf_swb[m,:,:]= dm9[i,:,:], dm9b[i,:,:]
  fluxtop_lwa[m,:,:], fluxtop_lwb[m,:,:]  = dm10[i,:,:], dm10b[i,:,:]
  fluxtop_swa[m,:,:], fluxtop_swb[m,:,:]  = dm11[i,:,:], dm11b[i,:,:]
  
  m = m+1

# Initialise variables to calculate differences between reference run and dust storm run
dustq_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
temp_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
ps_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
rho_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
u_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
 
fluxsurf_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxsurf_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))

# Calculate differences
for i in xrange(t.shape[0]):
 dustq_diff[i,:,:] = dustq_a[i,:,:] - dustq_b[i,:,:]
 temp_diff[i,:,:] = temp_a[i,:,:] - temp_b[i,:,:]
 ps_diff[i,:,:] = ps_a[i,:,:] - ps_b[i,:,:]
 rho_diff[i,:,:] = rho_a[i,:,:] - rho_b[i,:,:]
 u_diff[i,:,:] = u_a[i,:,:] - u_b[i,:,:]
  
 fluxsurf_lw_diff[i,:,:] = fluxsurf_lwa[i,:,:] - fluxsurf_lwb[i,:,:]
 fluxsurf_sw_diff[i,:,:] = fluxsurf_swa[i,:,:] - fluxsurf_swb[i,:,:]
 fluxtop_lw_diff[i,:,:] = fluxtop_lwa[i,:,:] - fluxtop_lwb[i,:,:]
 fluxtop_sw_diff[i,:,:] = fluxtop_swa[i,:,:] - fluxtop_swb[i,:,:]

##PLOTS

bwr = matplotlib.cm.bwr # Colormap that needs to be shifted set here

# Longitude slice
lon_slice = 8
print 'longitude selected: %.f' % lon[lon_slice]

# Zonal averaging
temp_avg = np.zeros((t.shape[0],lat.shape[0]))
for i in xrange(0,lon.shape[0],1):
 temp_avg = temp_avg + temp_diff[:,:,i]
 
temp_avg = temp_avg/len(lon)

#ps_t
#u_t
#v_t
#rho_t
#fssw_t
#fslw_t
#ftsw_t
#ftlw_t

# Time moving-point average

nn=200 # Number of points to average over

temp_avg_t = np.zeros((t.shape[0]-(nn-1),lat.shape[0]))

for i in xrange(0,lat.shape[0],1):
 temp_avg_t[:,i] = moving_average(temp_avg[:,i],n=nn)

# Temperature PLOT

temp_t = np.matrix.transpose(temp_avg_t) #np.matrix.transpose(temp_diff[:,:,lon_slice])

#print temp_t.shape,temp_avg.shape,temp_t,temp_avg

midd = (1 - np.amax(temp_t)/(np.amax(temp_t) + abs(np.amin(temp_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(temp_t,cmap=nbwr)
plt.axis('tight')
plt.title('Temperature difference (benchmark vs dust storm GCM runs) \n altitude: surface')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.f')
cb.set_label('Temperature / K', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/SurfTempDiff_LatvsTime_FY.png")

# Surface pressure PLOT

ps_t = np.matrix.transpose(ps_diff[:,:,lon_slice])

midd = (1 - np.amax(ps_t)/(np.amax(ps_t) + abs(np.amin(ps_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,ps_t,cmap=nbwr)
plt.axis('tight')
plt.title('Surface pressure difference (benchmark vs dust storm GCM runs) \n altitude: surface, Longitude: 140W')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.f')
cb.set_label('Surface pressure / Pa', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/SurfPresDiff_LatvsTime_FY.png") 

# Zonal wind PLOT

u_t = np.matrix.transpose(u_diff[:,:,lon_slice])

midd = (1 - np.amax(u_t)/(np.amax(u_t) + abs(np.amin(u_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,u_t,cmap=nbwr)
plt.axis('tight')
plt.title('Zonal wind speed difference (benchmark vs dust storm GCM runs) \n altitude: surface, Longitude: 140W')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.f')
cb.set_label('Zonal wind speed / ms^-1', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/ZonalWindDiff_LatvsTime_FY.png")

# Atmospheric density PLOT

rho_t = np.matrix.transpose(rho_diff[:,:,lon_slice])

midd = (1 - np.amax(rho_t)/(np.amax(rho_t) + abs(np.amin(rho_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,rho_t,cmap=nbwr)
plt.axis('tight')
plt.title('Atmospheric density difference (benchmark vs dust storm GCM runs) \n altitude: surface, Longitude: 140W')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1e')
cb.set_label('Atmospheric density / kg m^-3', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/DensityDiff_LatvsTime_FY.png")

# Longwave incoming radiation PLOT

fslw_t = np.matrix.transpose(fluxsurf_lw_diff[:,:,lon_slice])

midd = (1 - np.amax(fslw_t)/(np.amax(fslw_t) + abs(np.amin(fslw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,fslw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Incident long-wave (IR) surface flux (difference) \n altitude: surface, Longitude: 140W')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FSLW_Diff_LatvsTime_FY.png")

# Longwave outgoing radiation PLOT

ftlw_t = np.matrix.transpose(fluxtop_lw_diff[:,:,lon_slice])

midd = (1 - np.amax(ftlw_t)/(np.amax(ftlw_t) + abs(np.amin(ftlw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,ftlw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Outgoing long-wave (IR) surface flux (difference) \n altitude: surface, Longitude: 140W')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FTLW_Diff_LatvsTime_FY.png")

# Shortwave incident radiation PLOT

fssw_t = np.matrix.transpose(fluxsurf_sw_diff[:,:,lon_slice])

midd = (1 - np.amax(fssw_t)/(np.amax(fssw_t) + abs(np.amin(fssw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,fssw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Incident short-wave (IR) surface flux (difference) \n altitude: surface, Longitude: 140W')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FSSW_Diff_LatvsTime_FY.png")

# Shortwave outgoing radiation PLOT

ftsw_t = np.matrix.transpose(fluxtop_sw_diff[:,:,lon_slice])

midd = (1 - np.amax(ftsw_t)/(np.amax(ftsw_t) + abs(np.amin(ftsw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,ftsw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Outgoing short-wave (IR) surface flux (difference) \n altitude: surface, Longitude: 140W')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FTSW_Diff_LatvsTime_FY.png")
