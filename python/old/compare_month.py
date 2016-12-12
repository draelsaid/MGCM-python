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
from plt_timeseries import plt_timeseries, MidpointNormalize

rundira = "a_ds"
montha = "m5"
filenamea = "diagfi.nc"

rundirb = "a_ref"
monthb = "m5"
filenameb = "diagfi.nc"

a = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/MGCM_v5-1/runs/%s/%s/%s" % (rundira,montha,filenamea),'r')

b = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/MGCM_v5-1/runs/%s/%s/%s" % (rundirb,monthb,filenameb),'r')

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

c = np.matrix('1 45; 103 45; 243 2.5') # Dust storm points

major_ticksx = np.arange(0, 361, 40)                                              
minor_ticksx = np.arange(0, 361, 10)
major_ticksy = np.arange(-90, 91, 30)                                              
minor_ticksy = np.arange(-90, 91, 10)

## Temperature PLOT
temp_t = np.matrix.transpose(temp_diff[:,:,8])

fig = plt.figure(figsize=(10,10), dpi=100)
ax = fig.add_subplot(1,1,1)

plt.pcolormesh(t,lat,temp_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.xlabel('Solar longitude / degrees',fontsize=10)
plt.ylabel('Latitude / degrees',fontsize=10)

# Extra Markers
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 ax.plot(c[i,0],c[i,1],'o',color='y')

# Ticks
ax.set_xticks(major_ticksx)                                                       
ax.set_xticks(minor_ticksx, minor=True)                                       
ax.set_yticks(major_ticksy)                                                       
ax.set_yticks(minor_ticksy, minor=True)
ax.tick_params(axis='both', labelsize=10)

# Color bar
cb = plt.colorbar(format='%.0f')
cb.set_label('Temperature / K', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/SurfTempDiff_LatvsTime_FY_uavg_tavg.png")

## Surface pressure and Atmospheric density at surface PLOT
ps_t = np.matrix.transpose(ps_diff[:,:,8])
rho_t = np.matrix.transpose(rho_diff[:,:,8])

f, axarr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t
y = lat
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'

cb_label = 'Pressure / Pa'
cb_label2 = 'Density / kg / m^-3'
 
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')

ax1 = axarr[0].pcolormesh(x, y, ps_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
axarr[0].axis('tight') 
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 axarr[0].plot(c[i,0],c[i,1],'o',color='y')
axarr[0].set_xticks(major_ticksx)
axarr[0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0].set_yticks(major_ticksy)                                                       
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('(a)', fontsize=18)
axarr[0].tick_params(axis='both', labelsize=14)

ax2 = axarr[1].pcolormesh(x, y, rho_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 axarr[1].plot(c[i,0],c[i,1],'o',color='y')
axarr[1].set_title('(b)', fontsize=18)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.54, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax1, cax=cbar_ax, format='%.0f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=14)                    # colorbar label

cbar_ax2 = f.add_axes([0.85, 0.1, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb2 = f.colorbar(ax2, cax=cbar_ax2, format='%.1e', extend='both') # double-edged colorbar
cb2.set_label('%s' % (cb_label2), fontsize=14)                    # colorbar label

plt.savefig("/home/physastro/aes442/SurfPresDiff_SurfDensDiff_LatvsLs_zonavg_tavg.png") 
exit()





# Temperature PLOT

temp_t = np.matrix.transpose(temp_diff[:,:,8])

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,temp_t)
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
# Surface pressure and Atmospheric density at surface PLOT
ps_t = np.matrix.transpose(ps_diff[:,:,8])
rho_t = np.matrix.transpose(rho_diff[:,:,8])

f, axarr = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t
y = lat
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'

cb_label = 'Pressure / Pa'
cb_label2 = 'Density / kg / m^-3'
 
major_ticksx = np.arange(0, 361, 40)                                              
minor_ticksx = np.arange(0, 361, 10)
major_ticksy = np.arange(-90, 92, 30)                                              
minor_ticksy = np.arange(-90, 92, 10)
 
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')
   
ax1 = plt.pcolormesh(x, y, ps_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 f.gca().add_artist(cmarker)
 axarr[0,0].plot(c[i,0],c[i,1],'o',color='y')
axarr[0,0].set_xticks(major_ticksx)
axarr[0,0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0,0].set_yticks(major_ticksy)                                                       
axarr[0,0].set_yticks(minor_ticksy, minor=True)
axarr[0,0].set_title('(a)', fontsize=10)

ax2 = axarr[1,0].pcolormesh(x, y, rho_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 f.gca().add_artist(cmarker)
 axarr[1,0].plot(c[i,0],c[i,1],'o',color='y')
axarr[1,0].set_xticks(major_ticksx)                                                       
axarr[1,0].set_xticks(minor_ticksx, minor=True)                                           
axarr[1,0].set_yticks(major_ticksy)                                                       
axarr[1,0].set_yticks(minor_ticksy, minor=True)
axarr[1,0].set_title('(b)', fontsize=10)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.1, 0.05, 0.4])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax1, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=14)                    # colorbar label

cbar_ax2 = f.add_axes([0.85, 0.1, 0.05, 0.4])  # [h_placement, v_placement, h_size, v_size]
cb2 = f.colorbar(ax2, cax=cbar_ax2, format='%.1f', extend='both') # double-edged colorbar
cb2.set_label('%s' % (cb_label2), fontsize=14)                    # colorbar label

plt.savefig("/home/physastro/aes442/SurfPresDiff_SurfDensDiff_LatvsLs_zonavg_tavg.png") 

# Atmospheric density PLOT

rho_t = np.matrix.transpose(rho_diff[:,:,8])

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,rho_t)
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

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,fslw_t)
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

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,ftlw_t)
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

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,fssw_t)
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

plt.figure(figsize=(10,10))
plt.pcolormesh(t,lat,ftsw_t)
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
 
