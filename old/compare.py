# Compares NetCDF data from the Mars GCM
# Moving time graph code also included
# Adam El-Said 07/2016

#Time moving graph code
# plt.show(block=False)
# duration = 0.5
# plt.pause(duration)
# plt.close()

from scipy.io import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm,ticker
from matplotlib.ticker import LinearLocator

import matplotlib.pyplot as plt
import numpy as np
import time
import os,sys

from shiftedColorMap import *

rundira = "a_ds_writediag12"
montha = "m1"
filenamea = "diagfi.nc"

rundirb = "a_benchmark_writediag12"
monthb = "m1"
filenameb = "diagfi.nc"

a = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/runs/%s/%s/%s" % (rundira,montha,filenamea),'r')

b = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/runs/%s/%s/%s" % (rundirb,monthb,filenameb),'r')

lat = a.variables['lat'][:]
lon = a.variables['lon'][:]
sigma = a.variables['sigma'][:]
t = a.variables['time'][:]

psa = a.variables['ps'][:]
tempa = a.variables['temp'][:]
ua = a.variables['u'][:]
va = a.variables['v'][:]
dustqa = a.variables['dustq'][:]
dustNa = a.variables['dustN'][:]
rhoa = a.variables['rho'][:]
fluxsurf_lwa = a.variables['fluxsurf_lw'][:]
fluxsurf_swa = a.variables['fluxsurf_sw'][:]
fluxtop_lwa = a.variables['fluxtop_lw'][:]
fluxtop_swa = a.variables['fluxtop_sw'][:]

psb = b.variables['ps'][:]
tempb = b.variables['temp'][:]
ub = b.variables['u'][:]
vb = b.variables['v'][:]
dustqb = b.variables['dustq'][:]
dustNb = b.variables['dustN'][:]
rhob = b.variables['rho'][:]
fluxsurf_lwb = b.variables['fluxsurf_lw'][:]
fluxsurf_swb = b.variables['fluxsurf_sw'][:]
fluxtop_lwb = b.variables['fluxtop_lw'][:]
fluxtop_swb = b.variables['fluxtop_sw'][:]

dustq_a = dustqa[:,1,:,:]
dustq_b = dustqb[:,1,:,:]

temp_a = tempa[:,1,:,:]
temp_b = tempb[:,1,:,:]

rho_a = rhoa[:,1,:,:]
rho_b = rhob[:,1,:,:]

u_a = ua[:,1,:,:]
u_b = ub[:,1,:,:]

dustq_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
temp_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
ps_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
rho_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
u_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
 
fluxsurf_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxsurf_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
 
for i in xrange(t.shape[0]):
 dustq_diff[i,:,:] = dustq_a[i,:,:] - dustq_b[i,:,:]
 temp_diff[i,:,:] = temp_a[i,:,:] - temp_b[i,:,:]
 ps_diff[i,:,:] = psa[i,:,:] - psb[i,:,:]
 rho_diff[i,:,:] = rho_a[i,:,:] - rho_b[i,:,:]
 u_diff[i,:,:] = u_a[i,:,:] - u_b[i,:,:]
  
 fluxsurf_lw_diff[i,:,:] = fluxsurf_lwa[i,:,:] - fluxsurf_lwb[i,:,:]
 fluxsurf_sw_diff[i,:,:] = fluxsurf_swa[i,:,:] - fluxsurf_swb[i,:,:]
 fluxtop_lw_diff[i,:,:] = fluxtop_lwa[i,:,:] - fluxtop_lwb[i,:,:]
 fluxtop_sw_diff[i,:,:] = fluxtop_swa[i,:,:] - fluxtop_swb[i,:,:]
 
# Time Moving PLOT
#plt.figure(figsize=(10,10))
#tint = 4 #Time interval

#for i in xrange(t.shape[0]/tint):
# plt.contourf(lon,lat,temp_diff[i*tint,:,:])
# tsol = i*tint/float(4)
# plt.title('Temperature difference (benchmark vs dust storm GCM runs) \n altitude: surface \n sol: %.2f' % (tsol))#,y=1.02,x=0.65)
# plt.xlabel('Longtitude / degrees')
# plt.ylabel('Latitude / degrees')
 
# cb=plt.colorbar(format='%.2f',label='Temperature / K')
# tick_locator = ticker.MaxNLocator(nbins=16)
# cb.locator = tick_locator
# cb.update_ticks()

# plt.show()

bwr = matplotlib.cm.bwr # Colormap that needs to be shifted set here

# Temperature PLOT
temp_t = np.matrix.transpose(temp_diff[:,:,8])

midd = (1 - np.amax(temp_t)/(np.amax(temp_t) + abs(np.amin(temp_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,temp_t,cmap=nbwr)
plt.axis('tight')
plt.title('Temperature difference (benchmark vs dust storm GCM runs) \n altitude: surface')#,y=1.02,x=0.65)
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.f')
cb.set_label('Temperature / K', rotation=270)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/SurfTempDiff_LatvsTime_%s.png" % (montha))

# Surface pressure PLOT

ps_t = np.matrix.transpose(ps_diff[:,:,8])

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

plt.savefig("/home/physastro/aes442/SurfPresDiff_LatvsTime_%s.png" % (montha)) 

# Zonal wind PLOT

u_t = np.matrix.transpose(u_diff[:,:,8])

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

plt.savefig("/home/physastro/aes442/ZonalWindDiff_LatvsTime_%s.png" % (montha))

# Atmospheric density PLOT

rho_t = np.matrix.transpose(rho_diff[:,:,8])

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

plt.savefig("/home/physastro/aes442/DensityDiff_LatvsTime_%s.png" % (montha))

# Longwave incoming radiation PLOT

fslw_t = np.matrix.transpose(fluxsurf_lw_diff[:,:,8])

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

plt.savefig("/home/physastro/aes442/FSLW_Diff_LatvsTime_%s.png" % (montha))

# Longwave outgoing radiation PLOT

ftlw_t = np.matrix.transpose(fluxtop_lw_diff[:,:,8])

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

plt.savefig("/home/physastro/aes442/FTLW_Diff_LatvsTime_%s.png" % (montha))

# Shortwave incident radiation PLOT

fssw_t = np.matrix.transpose(fluxsurf_sw_diff[:,:,8])

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

plt.savefig("/home/physastro/aes442/FSSW_Diff_LatvsTime_%s.png" % (montha))

# Shortwave outgoing radiation PLOT
ftsw_t = np.matrix.transpose(fluxtop_sw_diff[:,:,8])

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

plt.savefig("/home/physastro/aes442/FTSW_Diff_LatvsTime_%s.png" % (montha))

 #plt.show(block=False)
 #duration = 0.5
 #plt.pause(duration)
 #plt.close()
 
