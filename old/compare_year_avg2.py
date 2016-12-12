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
fslwa, fslwb = {}, {}
fsswa, fsswb = {}, {}
ftlwa, ftlwb = {}, {}
ftswa, ftswb = {}, {}

t = np.arange(1,669,float(1)/4) # Create new time dimension variable. Be careful to change interval if writediagfi interval is changed in runscript
 
# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,13):
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
 
 psa[i], psb[i] = a.variables['ps'][:], b.variables['ps'][:]
 tempa[i], tempb[i] = a.variables['temp'][:], b.variables['temp'][:]
 ua[i], ub[i] = a.variables['u'][:], b.variables['u'][:]
 dustqa[i], dustqb[i] = a.variables['dustq'][:], b.variables['dustq'][:]
 dustNa[i], dustNb[i] = a.variables['dustN'][:], b.variables['dustN'][:]
 rhoa[i], rhob[i] = a.variables['rho'][:], b.variables['rho'][:]
 fslwa[i], fslwb[i] = a.variables['fluxsurf_lw'][:], b.variables['fluxsurf_lw'][:]
 fsswa[i], fsswb[i] = a.variables['fluxsurf_sw'][:], b.variables['fluxsurf_sw'][:]
 ftlwa[i], ftlwb[i] = a.variables['fluxtop_lw'][:], b.variables['fluxtop_lw'][:]
 ftswa[i], ftswb[i] = a.variables['fluxtop_sw'][:], b.variables['fluxtop_sw'][:]

# Get time dimension length
n = 0
for i in xrange(1,len(psa)+1,1): # len(psa) gives the number of months
 n = n + len(dustqa[i])          # len(dustqa[i]) gives the number of time steps in each month. Different variable used as a cross-check of dimension consistency.

# dictionary.get(x,y) gives you the xth entry in the dictionary
# dictionary.keys() givs you the array references each dict entry

# Calculate differences
ps_d = {key: psa[key] - psb.get(key, 0) for key in psa.keys()}
dustq_d = {key: dustqa[key] - dustqb.get(key, 0) for key in dustqa.keys()}
temp_d = {key: tempa[key] - tempb.get(key, 0) for key in tempa.keys()}
rho_d = {key: rhoa[key] - rhob.get(key, 0) for key in rhoa.keys()}
u_d = {key: ua[key] - ub.get(key, 0) for key in ua.keys()}
fslw_d = {key: fslwa[key] - fslwb.get(key, 0) for key in fslwa.keys()}
fssw_d = {key: fsswa[key] - fsswb.get(key, 0) for key in fsswa.keys()}
ftlw_d = {key: ftlwa[key] - ftlwb.get(key, 0) for key in ftlwa.keys()}
ftsw_d = {key: ftswa[key] - ftswb.get(key, 0) for key in ftswa.keys()}

print ps_d
exit()

# Create numpy arrays out of the dictionaries
ps_a, ps_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
temp_a, temp_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
u_a, u_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustq_a, dustq_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustN_a, dustN_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
rho_a, rho_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_lwa, fluxsurf_lwb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_swa, fluxsurf_swb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_lwa, fluxtop_lwb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_swa, fluxtop_swb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))

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
