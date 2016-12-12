# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import math

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
Ls_m = {}
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

t = np.arange(1,669,float(1)/12) # Create new time dimension variable. Be careful to change interval if writediagfi interval is changed in runscript
 
# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,13,1):
 rundira = "a_ds_writediag12"
 rundirb = "a_benchmark_writediag12" 
 month = "m%i" % (i)
 filename = "diagfi.nc"
 
 a = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/runs/%s/%s/%s" % (rundira,month,filename),'r')
 b = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/runs/%s/%s/%s" % (rundirb,month,filename),'r')

 lat = a.variables['lat'][:]
 lon = a.variables['lon'][:]
 sigma = a.variables['sigma'][:]
 t_m = a.variables['time'][:]
 Ls_m[i] = a.variables['Ls'][:]
 
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

# Create solar longitude vector (so it fits with sol vector - bit scruffy)
Ls = np.zeros((n))

p=0
for i in xrange(1,len(Ls_m)+1,1):
 gg = Ls_m[i]
 for j in xrange(gg.shape[0]):
  Ls[p] = gg[j]
  p = p + 1

Ls = Ls[:-12]

exit()

# Create all other variables, with altitude dimension removed
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

# final loop to allocate each run to a variable (numpy arrays) 
 for i in xrange(0,dm2.shape[0],1):
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
  
  m = m+1

tempaa, tempbb = {}, {}
uaa, ubb = {}, {}
dustqaa, dustqbb = {}, {}
dustNaa, dustNbb = {}, {}
rhoaa, rhobb = {}, {}

for k in xrange(1,len(psa)+1,1):
 dm2, dm2b = tempa[k], tempb[k]
 dm3, dm3b = ua[k], ub[k]
 dm5, dm5b = dustqa[k], dustqb[k]
 dm6, dm6b = dustNa[k], dustNb[k]
 dm7, dm7b = rhoa[k], rhob[k]

 d2, d2b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0]))
 d3, d3b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0]))
 d5, d5b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0]))
 d6, d6b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0]))
 d7, d7b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0]))

 for j in xrange(0,lon.shape[0],1):
  d2, d2b = d2 + dm2[:,:,:,j], d2b + dm2b[:,:,:,j]
  d3, d3b = d3 + dm3[:,:,:,j], d3b + dm3b[:,:,:,j]
  d5, d5b = d5 + dm5[:,:,:,j], d5b + dm5b[:,:,:,j]
  d6, d6b = d6 + dm6[:,:,:,j], d6b + dm6b[:,:,:,j]
  d7, d7b = d7 + dm7[:,:,:,j], d7b + dm7b[:,:,:,j]
  
 tempaa[k], tempbb[k] = d2, d2b
 uaa[k], ubb[k] = d3, d3b
 dustqaa[k], dustqbb[k] = d5, d5b
 dustNaa[k], dustNbb[k] = d6, d6b
 rhoaa[k], rhobb[k] = d7, d7b

# Create variables now without longitude
temp_aa, temp_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
u_aa, u_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustq_aa, dustq_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustN_aa, dustN_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rho_aa, rho_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))

m=0
for k in xrange(1,len(psa)+1,1):
 dd2, dd2b = tempaa[k], tempbb[k]
 dd3, dd3b = uaa[k], ubb[k]
 dd5, dd5b = dustqaa[k], dustqbb[k]
 dd6, dd6b = dustNaa[k], dustNbb[k]
 dd7, dd7b = rhoaa[k], rhobb[k]

 for i in xrange(dd2.shape[0]):
  temp_aa[m,:,:], temp_bb[m,:,:]   = dd2[i,:,:], dd2b[i,:,:]
  u_aa[m,:,:], u_bb[m,:,:]         = dd3[i,:,:], dd3b[i,:,:]
  dustq_aa[m,:,:], dustq_bb[m,:,:] = dd5[i,:,:], dd5b[i,:,:]
  dustN_aa[m,:,:], dustN_bb[m,:,:] = dd6[i,:,:], dd6b[i,:,:]
  rho_aa[m,:,:], rho_bb[m,:,:]     = dd7[i,:,:], dd7b[i,:,:]

  m = m+1

# Finally create the zonal averaged data (time,alt,lat)
temp_aa, temp_bb   = temp_aa/len(lon), temp_bb/len(lon)
u_aa, u_bb         = u_aa/len(lon), u_bb/len(lon) 
dustq_aa, dustq_bb = dustq_aa/len(lon), dustq_bb/len(lon)
dustN_aa, dustN_bb = dustN_aa/len(lon), dustN_bb/len(lon) 
rho_aa, rho_bb     = rho_aa/len(lon), rho_bb/len(lon) 

# Initialise variables to calculate differences between reference run and dust storm run
dustq_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
dustN_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
temp_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
ps_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
rho_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
u_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
 
fluxsurf_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxsurf_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))

dq_d = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
dN_d = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
t_d = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
rho_d = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
u_d = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))

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
 
 t_d[i,:,:] = temp_aa[i,:,:] - temp_bb[i,:,:]
 dq_d[i,:,:] = dustq_aa[i,:,:] - dustq_bb[i,:,:]
 dN_d[i,:,:] = dustN_aa[i,:,:] - dustN_bb[i,:,:]
 rho_d[i,:,:] = rho_aa[i,:,:] - rho_bb[i,:,:]
 u_d[i,:,:] = u_aa[i,:,:] - u_bb[i,:,:]

# Zonal averaging (time,lat)
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

 ########### NOW NEED TIME AVERAGE FOR NEW DATA ################### (SEASONAL)
nnn=240 # Number of points to average over
t_av = t[:-(nnn-1)]

td_avg = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
dqd_avg = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
dNd_avg = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
rhod_avg = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
ud_avg = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))

for j in xrange(0,lat.shape[0],1):
 for i in xrange(0,sigma.shape[0],1):
  td_avg[:,i,j] = moving_average(t_d[:,i,j],n=nnn)
  dqd_avg[:,i,j] = moving_average(dq_d[:,i,j],n=nnn)
  dNd_avg[:,i,j] = moving_average(dN_d[:,i,j],n=nnn)
  rhod_avg[:,i,j] = moving_average(rho_d[:,i,j],n=nnn)
  ud_avg[:,i,j] = moving_average(u_d[:,i,j],n=nnn)

# Calculate approximate HEIGHT from sigma (km)
alt = np.zeros((sigma.shape[0]))
for i in xrange(len(sigma)):
 alt[i] = -10.8*math.log(sigma[i])

# Temperature PLOT
temp_t = np.matrix.transpose(temp_avg_t)

plt.figure()
plt.pcolormesh(t_avg,lat,temp_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
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

plt.figure()
plt.pcolormesh(t_avg,lat,ps_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Surface pressure / Pa', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/SurfPresDiff_LatvsTime_FY_uavg_tavg.png") 

# Zonal wind PLOT
u_t = np.matrix.transpose(u_avg_t)

plt.figure()
plt.pcolormesh(t_avg,lat,u_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Zonal wind velocity / ms^-1', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/ZonalWindDiff_LatvsTime_FY_uavg_tavg.png")

# Atmospheric density PLOT
rho_t = np.matrix.transpose(rho_avg_t)

plt.figure()
plt.pcolormesh(t_avg,lat,rho_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1e')
cb.set_label('Atmospheric density / kg m^-3', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/DensityDiff_LatvsTime_FY_uavg_tavg.png")

# Longwave incoming radiation PLOT
fslw_t = np.matrix.transpose(fslw_avg_t)

plt.figure()
plt.pcolormesh(t_avg,lat,fslw_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FSLW_Diff_LatvsTime_FY_uavg_tavg.png")

# Longwave outgoing radiation PLOT
ftlw_t = np.matrix.transpose(ftlw_avg_t)

plt.figure()
plt.pcolormesh(t_avg,lat,ftlw_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FTLW_Diff_LatvsTime_FY_uavg_tavg.png")

# Shortwave incident radiation PLOT
fssw_t = np.matrix.transpose(fssw_avg_t)

plt.figure()
plt.pcolormesh(t_avg,lat,fssw_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FSSW_Diff_LatvsTime_FY_uavg_tavg.png")

# Shortwave outgoing radiation PLOT
ftsw_t = np.matrix.transpose(ftsw_avg_t)

plt.figure()
plt.pcolormesh(t_avg,lat,ftsw_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Time / sols')
plt.ylabel('Latitude / degrees')
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Surface flux / W m^-2', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/FTSW_Diff_LatvsTime_FY_uavg_tavg.png")

## Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff, t, lat, lon, 5, 4, 5, 'Longitude / degrees', 'Latitude / degrees', 'Sol', 'Dust MMR / kg/kg', 36, 'Dustq_zontimeavg_tseries.png')

### Exo Mars
## Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff[5628:,:,:], t[5628:], lat, lon, 5, 4, 5, 'Longitude / degrees', 'Latitude / degrees', 'Sol', 'Dust MMR / kg/kg', 36, 'Dustq_zontimeavg_tseries_exomars.png')

## Time series temperature (K) (time, alt, lat)
plt_timeseries(td_avg[:,1:31,:], t_av, alt[1:31], lat, 2, 2, 5, 'Altitude / km', 'Latitude / degrees', 'Season', 'Temperature / K', 1940, 'TempDiff_zontimeavg_tseries.png')

## Time series atmospheric density (time, alt, lat)
plt_timeseries(rhod_avg[:,1:31,:], t_av, alt[1:31], lat, 2, 2, 5, 'Altitude / km', 'Latitude / degrees', 'Season', 'Atmospheric density / kg m^-3', 1940, 'Density_zontimeavg_tseries.png')

## Time series zonal wind velocity (time, alt, lat)
plt_timeseries(ud_avg[:,1:31,:], t_av, alt[1:31], lat, 2, 2, 5, 'Altitude / km', 'Latitude / degrees', 'Season', 'Zonal wind velocity / ms^-1', 1940, 'ZonalWindVel_zontimeavg_tseries.png')

#fix trigger to change from contour to pcolormesh
#fix trigger to label plots properly
#fix trigger to change the magnitude of units in colorbar


