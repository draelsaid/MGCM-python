# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

from shiftedColorMap import *
from scipy.io import *
from matplotlib import cm,ticker
from plt_timeseries import plt_timeseries, MidpointNormalize

# Moving average 
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Initialise dictionaries - due to data size
psa, psb = {}, {} 
tempa, tempb = {}, {} 
ua, ub = {}, {}
dustqa, dustqb = {}, {}
dustNa, dustNb = {}, {}
rhoa, rhob = {}, {}
fluxsurflwa, fluxsurflwb = {}, {}
fluxsurfswa, fluxsurfswb = {}, {}
fluxtoplwa, fluxtoplwb = {}, {}
fluxtopswa, fluxtopswb = {}, {}

t = np.arange(1,669,float(1)/4) # Create new time dimension variable. Be careful to change interval if writediagfi interval is changed in runscript
 
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
 #Ls = a.variables['Ls'][:]
 
 psa[i] = a.variables['ps'][:]
 tempa[i] = a.variables['temp'][:]
 ua[i] = a.variables['u'][:]
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

# Create all other variables, with altitude dimension removed, to save time
ps_a, ps_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
temp_a, temp_b   = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
u_a, u_b         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustq_a, dustq_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustN_a, dustN_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
rho_a, rho_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_lwa, fluxsurf_lwb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_swa, fluxsurf_swb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_lwa, fluxtop_lwb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_swa, fluxtop_swb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))

# For zonally averaged (time,alt,lat) data
temp_aa, temp_bb   = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
u_aa, u_bb         = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustq_aa, dustq_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustN_aa, dustN_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rho_aa, rho_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))

d2, d2b = np.zeros((244,sigma.shape[0],lat.shape[0])), np.zeros((244,sigma.shape[0],lat.shape[0]))
d3, d3b = np.zeros((244,sigma.shape[0],lat.shape[0])), np.zeros((244,sigma.shape[0],lat.shape[0]))
d5, d5b = np.zeros((244,sigma.shape[0],lat.shape[0])), np.zeros((244,sigma.shape[0],lat.shape[0]))
d6, d6b = np.zeros((244,sigma.shape[0],lat.shape[0])), np.zeros((244,sigma.shape[0],lat.shape[0]))
d7, d7b = np.zeros((244,sigma.shape[0],lat.shape[0])), np.zeros((244,sigma.shape[0],lat.shape[0]))

m=0
for k in xrange(1,len(psa)+1,1):
# Dummy variables - extract month k from dictionary
 dm1, dm1b = psa[k], psb[k]
 dm2, dm2b = tempa[k], tempb[k]
 dm3, dm3b = ua[k], ub[k]
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
 dmm5, dmm5b = dm5[:,1,:,:], dm5b[:,1,:,:]
 dmm6, dmm6b = dm6[:,1,:,:], dm6b[:,1,:,:]
 dmm7, dmm7b = dm7[:,1,:,:], dm7b[:,1,:,:] 
 
# Zonally averaged (time,sigma,lat) data
 for k in xrange(len(lon)):
  d2, d2b = d2 + dm2[:,:,:,k], d2 + dm2b[:,:,:,k]
  d3, d3b = d3 + dm3[:,:,:,k], d3 + dm3b[:,:,:,k]
  d5, d5b = d5 + dm5[:,:,:,k], d5 + dm5b[:,:,:,k]
  d6, d6b = d6 + dm6[:,:,:,k], d6 + dm6b[:,:,:,k]
  d7, d7b = d7 + dm7[:,:,:,k], d7 + dm7b[:,:,:,k]
  
 d2, d2b = d2/len(lon), d2b/len(lon)
 d3, d3b = d3/len(lon), d3b/len(lon)
 d5, d5b = d5/len(lon), d5b/len(lon)
 d6, d6b = d6/len(lon), d6b/len(lon)
 d7, d7b = d7/len(lon), d7b/len(lon)

# final loop to allocate each run to a variable (numpy arrays) 
 for i in xrange(dm2.shape[0]):
  ps_a[m,:,:], ps_b[m,:,:]                = dm1[i,:,:], dm1b[i,:,:]
  temp_a[m,:,:], temp_b[m,:,:]            = dmm2[i,:,:], dmm2b[i,:,:]
  u_a[m,:,:], u_b[m,:,:]                  = dmm3[i,:,:], dmm3b[i,:,:]
  dustq_a[m,:,:], dustq_b[m,:,:]          = dmm5[i,:,:], dmm5b[i,:,:]
  dustN_a[m,:,:], dustN_b[m,:,:]          = dmm6[i,:,:], dmm6b[i,:,:]
  rho_a[m,:,:], rho_b[m,:,:]              = dmm7[i,:,:], dmm7b[i,:,:]
  fluxsurf_lwa[m,:,:], fluxsurf_lwb[m,:,:]= dm8[i,:,:], dm8b[i,:,:]
  fluxsurf_swa[m,:,:], fluxsurf_swb[m,:,:]= dm9[i,:,:], dm9b[i,:,:]
  fluxtop_lwa[m,:,:], fluxtop_lwb[m,:,:]  = dm10[i,:,:], dm10b[i,:,:]
  fluxtop_swa[m,:,:], fluxtop_swb[m,:,:]  = dm11[i,:,:], dm11b[i,:,:]
  
  temp_aa[m,:,:], temp_bb[m,:,:]            = d2[i,:,:], d2b[i,:,:]
  u_aa[m,:,:], u_bb[m,:,:]                  = d3[i,:,:], d3b[i,:,:]
  dustq_aa[m,:,:], dustq_bb[m,:,:]          = d5[i,:,:], d5b[i,:,:]
  dustN_aa[m,:,:], dustN_bb[m,:,:]          = d6[i,:,:], d6b[i,:,:]
  rho_aa[m,:,:], rho_bb[m,:,:]              = d7[i,:,:], d7b[i,:,:]
  
  m = m+1

# Initialise variables to calculate differences between reference run and dust storm run
dustq_diff, dustN_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0])), np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
temp_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
ps_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
rho_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
u_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
 
fluxsurf_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxsurf_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))

dq_d = np.zeros((n,sigma.shape[0],lat.shape[0]))
dN_d = np.zeros((n,sigma.shape[0],lat.shape[0]))
t_d = np.zeros((n,sigma.shape[0],lat.shape[0]))
rho_d = np.zeros((n,sigma.shape[0],lat.shape[0]))
u_d = np.zeros((n,sigma.shape[0],lat.shape[0]))

# Calculate differences
for i in xrange(t.shape[0]):
 dustq_diff[i,:,:] = dustq_a[i,:,:] - dustq_b[i,:,:]
 dustN_diff[i,:,:] = dustN_a[i,:,:] - dustN_b[i,:,:]
 temp_diff[i,:,:] = temp_a[i,:,:] - temp_b[i,:,:]
 ps_diff[i,:,:] = ps_a[i,:,:] - ps_b[i,:,:]
 rho_diff[i,:,:] = rho_a[i,:,:] - rho_b[i,:,:]
 u_diff[i,:,:] = u_a[i,:,:] - u_b[i,:,:]
  
 fluxsurf_lw_diff[i,:,:] = fluxsurf_lwa[i,:,:] - fluxsurf_lwb[i,:,:]
 fluxsurf_sw_diff[i,:,:] = fluxsurf_swa[i,:,:] - fluxsurf_swb[i,:,:]
 fluxtop_lw_diff[i,:,:] = fluxtop_lwa[i,:,:] - fluxtop_lwb[i,:,:]
 fluxtop_sw_diff[i,:,:] = fluxtop_swa[i,:,:] - fluxtop_swb[i,:,:]
 
 # Zonally averaged (time,sigma,lat) data
 dq_d  = dustq_aa[i,:,:] - dustq_bb[i,:,:]
 dN_d  = dustN_aa[i,:,:] - dustN_bb[i,:,:] 
 t_d   = temp_aa[i,:,:] - temp_bb[i,:,:]
 rho_d = rho_aa[i,:,:] - rho_bb[i,:,:]
 u_d   = u_aa[i,:,:] - u_bb[i,:,:]

# Zonal averaging
temp_avg = np.zeros((t.shape[0],lat.shape[0]))
ps_avg = np.zeros((t.shape[0],lat.shape[0]))
u_avg = np.zeros((t.shape[0],lat.shape[0]))
rho_avg = np.zeros((t.shape[0],lat.shape[0]))
fssw_avg = np.zeros((t.shape[0],lat.shape[0]))
fslw_avg = np.zeros((t.shape[0],lat.shape[0]))
ftsw_avg = np.zeros((t.shape[0],lat.shape[0]))
ftlw_avg = np.zeros((t.shape[0],lat.shape[0]))

for i in xrange(0,lon.shape[0],1):
 temp_avg = temp_avg + temp_diff[:,:,i]
 ps_avg = ps_avg + ps_diff[:,:,i]
 u_avg = u_avg + u_diff[:,:,i]
 rho_avg = rho_avg + rho_diff[:,:,i]
 fssw_avg = fssw_avg + fluxsurf_sw_diff[:,:,i]
 fslw_avg = fslw_avg + fluxsurf_lw_diff[:,:,i]
 ftsw_avg = ftsw_avg + fluxtop_sw_diff[:,:,i]
 ftlw_avg = ftlw_avg + fluxtop_lw_diff[:,:,i]
 
temp_avg = temp_avg/len(lon)
ps_avg = ps_avg/len(lon)
u_avg = u_avg/len(lon)
rho_avg = rho_avg/len(lon)
fssw_avg = fssw_avg/len(lon)
fslw_avg = fslw_avg/len(lon)
ftsw_avg = ftsw_avg/len(lon)
ftlw_avg = ftlw_avg/len(lon)

# Time moving-point average OF Zonal average

nn=60 # Number of points to average over
t_avg = t[:-(nn-1)]

temp_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
ps_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
u_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
rho_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
fssw_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
fslw_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
ftsw_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
ftlw_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))

for i in xrange(0,lat.shape[0],1):
 temp_avg_t[:,i] = moving_average(temp_avg[:,i],n=nn)
 ps_avg_t[:,i] = moving_average(ps_avg[:,i],n=nn)
 u_avg_t[:,i] = moving_average(u_avg[:,i],n=nn)
 rho_avg_t[:,i] = moving_average(rho_avg[:,i],n=nn)
 fssw_avg_t[:,i] = moving_average(fssw_avg[:,i],n=nn)
 fslw_avg_t[:,i] = moving_average(fslw_avg[:,i],n=nn)
 ftsw_avg_t[:,i] = moving_average(ftsw_avg[:,i],n=nn)
 ftlw_avg_t[:,i] = moving_average(ftlw_avg[:,i],n=nn)
 
##PLOTS
bwr = matplotlib.cm.bwr # Colormap that needs to be shifted set here

# Longitude slice
lon_slice = 8
print 'longitude selected: %.f' % lon[lon_slice]
 
# Temperature PLOT
temp_t = np.matrix.transpose(temp_avg_t)

midd = (1 - np.amax(temp_t)/(np.amax(temp_t) + abs(np.amin(temp_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,temp_t,cmap=nbwr)
plt.axis('tight')
plt.title('Surface temperature difference \n zonal mean and moving time-average (%.i sols)' % (nn))
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Temperature / K', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/SurfTempDiff_LatvsTime_FY_uavg_tavg.png")

# Surface pressure PLOT
ps_t = np.matrix.transpose(ps_avg_t)

midd = (1 - np.amax(ps_t)/(np.amax(ps_t) + abs(np.amin(ps_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,ps_t,cmap=nbwr)
plt.axis('tight')
plt.title('Surface pressure difference \n zonal mean and moving time-average (%.i sols)')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Surface pressure / Pa', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/SurfPresDiff_LatvsTime_FY_uavg_tavg.png") 

# Zonal wind PLOT
u_t = np.matrix.transpose(u_avg_t) #np.matrix.transpose(u_diff[:,:,lon_slice])

midd = (1 - np.amax(u_t)/(np.amax(u_t) + abs(np.amin(u_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,u_t,cmap=nbwr)
plt.axis('tight')
plt.title('Surface zonal wind velocity difference \n zonal mean and moving time-average (%.i sols)')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Zonal wind speed / ms^-1', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/ZonalWindDiff_LatvsTime_FY_uavg_tavg.png")

# Atmospheric density PLOT

rho_t = np.matrix.transpose(rho_avg_t) #np.matrix.transpose(rho_diff[:,:,lon_slice])

midd = (1 - np.amax(rho_t)/(np.amax(rho_t) + abs(np.amin(rho_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,rho_t,cmap=nbwr)
plt.axis('tight')
plt.title('Surface atmospheric density difference \n zonal mean and moving time-average (%.i sols)')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1e')
cb.set_label('Atmospheric density / kg m^-3', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/DensityDiff_LatvsTime_FY_uavg_tavg.png")

# Longwave incoming radiation PLOT

fslw_t = np.matrix.transpose(fslw_avg_t) #np.matrix.transpose(fluxsurf_lw_diff[:,:,lon_slice])

midd = (1 - np.amax(fslw_t)/(np.amax(fslw_t) + abs(np.amin(fslw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,fslw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Incident LW (IR) surface flux difference \n zonal mean and moving time-average (%.i sols)')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FSLW_Diff_LatvsTime_FY_uavg_tavg.png")

# Longwave outgoing radiation PLOT

ftlw_t = np.matrix.transpose(ftlw_avg_t) #np.matrix.transpose(fluxtop_lw_diff[:,:,lon_slice])

midd = (1 - np.amax(ftlw_t)/(np.amax(ftlw_t) + abs(np.amin(ftlw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,ftlw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Outgoing LW (IR) surface flux difference \n zonal mean and moving time-average (%.i sols)')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FTLW_Diff_LatvsTime_FY_uavg_tavg.png")

# Shortwave incident radiation PLOT
fssw_t = np.matrix.transpose(fssw_avg_t)

midd = (1 - np.amax(fssw_t)/(np.amax(fssw_t) + abs(np.amin(fssw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,fssw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Incident SW (IR) surface flux, difference plot \n altitude: surface')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FSSW_Diff_LatvsTime_FY_uavg_tavg.png")

# Shortwave outgoing radiation PLOT
ftsw_t = np.matrix.transpose(ftsw_avg_t)

midd = (1 - np.amax(ftsw_t)/(np.amax(ftsw_t) + abs(np.amin(ftsw_t))))
nbwr = shiftedColorMap(bwr, midpoint=midd)

plt.figure()
plt.pcolormesh(t_avg,lat,ftsw_t,cmap=nbwr)
plt.axis('tight')
plt.title('Outgoing SW (IR) surface flux difference \n zonal mean and moving time-average (%.i sols)')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FTSW_Diff_LatvsTime_FY_uavg_tavg.png")

## Time series dust mmr
plt_timeseries(dustq_diff, t, lat, lon, 5, 4, 5, 'Longitude / degrees', 'Latitude / degrees', 'Sol', 'Dust MMR / kg/kg', 12, 'dustq_time_series.png')

## Time series dust mmr
plt_timeseries(temp_avg, t, lat, lon, 2, 2, 5, 'Latitude / degrees', 'Altitude / km', 'Sol', 'Dust MMR / kg/kg', 12, 'dustq_time_series.png')

## Want to now create LAT vs ALT temp_diff plots, zonally averaged, showing seasonal changes.
## Need to change the way the 4D data above is reduced to 3D (by selecting surface).
