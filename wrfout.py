# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016

import matplotlib as mpl
#mpl.use('Agg') # removes need for X-Server (graphics in linux). For qsub only.

import numpy as np
import pylab as py
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import datetime

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mars_time import MarsTime
from scipy.io import *
from matplotlib import cm,ticker
from plt_timeseries import *
from matplotlib.ticker import FormatStrFormatter
from MidPointNorm import *

# Prints EVERYTHING inside a variable without holding back (intended for diagnostic)
np.set_printoptions(threshold=np.inf)

# Abbreviate sol_ls conversion function
sol_Ls=MarsTime().sol_ls

# Moving average 
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# Grab topography from surface.nc or mola32.nc file
ml = netcdf.netcdf_file('/padata/alpha/users/aes442/datafile/surface.nc','r')

mola = {}

mola[0] = ml.variables['latitude'][:]
mola[1] = ml.variables['longitude'][:]
mola[2] = ml.variables['zMOL'][:]

# Import data from Luca's TES dust files for comparison
a = netcdf.netcdf_file('/padata/alpha/users/aes442/datafile/dust_MY28.nc','r')

d_lat_s = a.variables['latitude'][:]
d_lon_s = a.variables['longitude'][:]
d_t   = a.variables['Time'][:]
d_d   = a.variables['dustop'][:]

d_lat = np.linspace(-90,90,d_lat_s.shape[0])
d_lon = np.linspace(-180,180,d_lon_s.shape[0])

# Initialise dictionaries - due to data size
Ls_m = {}
psa, psb = {}, {} 
presa, presb = {}, {}
tempa, tempb = {}, {} 
tsurfa, tsurfb = {}, {} 
ua, ub = {}, {}
va, vb = {}, {}
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
dqsseda, dqssedb   = {}, {}
dqsdeva, dqsdevb   = {}, {}

# Number of files in run
f_files = np.array([6,7,8,9,10,11,12])

# DATA SHAPE >>> 
# bottom_top = ALTITUDE
# south_north= LATITUDE
# west_east  = LONGITUDE
#                   4D: (Time, bottom_top, south_north, west_east)
#                   3D: (Time, south_north, west_east)

for i in f_files:
 mgcm = "LMD_MMGCM"
 rundira = "a_ref"
 rundirb = "a_ds8"
 filename = "wrfout_d01_2024-01-0%i_00:00:00" % (i)
 
 a = netcdf.netcdf_file("/padata/alpha/users/aes442/mars_gcms/RUNS/R-%s/%s/%s/%s" % (mgcm,rundira,filename),'r')
 b = netcdf.netcdf_file("/padata/alpha/users/aes442/mars_gcms/RUNS/R-%s/%s/%s/%s" % (mgcm,rundirb,filename),'r')

## Referencing spatial/temporal variables
# Altitude
 znu      = a.variables['ZNU'][:]        
 znw      = a.variables['ZNW'][:]
 hgt[i]   = a.variables['HGT'][:]
 phtot[i] = a.variables['PHTOT'][:]

# Latitude and Longitude
 xlat = a.variables['XLAT'][:][0,:,0]
 xlon = a.variables['XLONG'][0,0,:]

# Required variables 
 dusta[i]  = a.variables['QDUST'][:]
 ptota[i]  = a.variables['PTOT'][:]
 tempa[i]  = a.variables['T'][:]
 ua[i]     = a.variables['U'][:]
 va[i]     = a.variables['V'][:]
 wa[i]     = a.variables['W'][:]
 taua[i]   = a.variables['TAU_DUST'][:]

 dustb[i]  = b.variables['QDUST'][:]
 ptotb[i]  = b.variables['PTOT'][:]
 tempb[i]  = b.variables['T'][:]
 ub[i]     = b.variables['U'][:]
 vb[i]     = b.variables['V'][:]
 wb[i]     = b.variables['W'][:]
 taub[i]   = b.variables['TAU_DUST'][:]

# Calculate approximate HEIGHT from geopotential
alt = np.zeros((sigma.shape[0]))

alt = phtot[0,:,0,0] / 3.72 - hgt[0,0,0]

print "Latitude: %i, (%.2f - %.2f) || Longitude: %i, (%.2f - %.2f) || Model levels: %i => Alt Min:%.3f | Alt Max:%.3f " % (xlat.shape[0],xlat[0],xlat[-1],xlon.shape[0],xlon[0],xlon[-1],alt.shape[0],alt[0],alt[-1])

# Time



fpath = "/home/physastro/aes442/results/Dustruns/m%i/" % (amth)

#########################################################################################

## PLOTS

# Common settings (ticks)
t_t = np.linspace(Ls_ee,Ls_e,n)
t_tau = np.linspace(Ls_ee,Ls_e,n)
lat_t = np.linspace(90,-90,lat.shape[0])
lon_t = np.linspace(-180,180,lon.shape[0])

# Solar longitude
i_mj=0.2
i_mn=0.05
major_ticksx = np.arange(Ls_ee, Ls_e+i_mj, i_mj)                                              
minor_ticksx = np.arange(Ls_ee, Ls_e, i_mn)
# Latitude
major_ticksy = np.arange(-90, 91, 30)                                              
minor_ticksy = np.arange(-90, 91, 10)

## tau_ref, tau_ds, tau_tes PLOT
tau_ds = np.matrix.transpose(tau_a_avg)
tau_ds =  tau_ds[:,l_1:l_2]
tau_ref = np.matrix.transpose(tau_b_avg)
tau_ref = tau_ref[:,l_1:l_2]
tau_TES = np.matrix.transpose(tau_d_avg)

f, axarr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t_tau
y = lat_t
xx = d_Ls_avg
yy = np.linspace(-90,90,d_lat_s.shape[0])
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'

cb_label = 'Dust optical depth / SI'
 
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')

ax1 = axarr[0].pcolormesh(x, y, tau_ds, cmap='gist_rainbow_r', vmin=np.min((np.min(tau_ds),np.min(tau_ref),np.min(tau_TES))), vmax=np.max((np.max(tau_ds),np.max(tau_ref),np.max(tau_TES))))
axarr[0].axis('tight') 
axarr[0].plot(c[0,0],c[0,1],'o',color='y',markersize=10)
axarr[0].set_xticks(major_ticksx)
axarr[0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0].set_yticks(major_ticksy)                                                       
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('(a) Dust storm run', fontsize=14)
axarr[0].tick_params(axis='both', labelsize=11, pad=10)

ax2 = axarr[1].pcolormesh(x, y, tau_ref, cmap='gist_rainbow_r', vmin=np.min((np.min(tau_ds),np.min(tau_ref),np.min(tau_TES))), vmax=np.max((np.max(tau_ds),np.max(tau_ref),np.max(tau_TES))))
axarr[1].set_title('(b) Reference run', fontsize=14)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8])                    # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax1, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=16)                    # colorbar label

#f.subplots_adjust(right=0.8)
#cbar_ax = f.add_axes([0.85, 0.665, 0.04, 0.235])  # [h_placement, v_placement, h_size, v_size]
#cb = f.colorbar(ax1, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
#cb.set_label('%s' % (cb_label), fontsize=16)                    # colorbar label

#f.subplots_adjust(right=0.8)
#cbar_ax2 = f.add_axes([0.85, 0.38, 0.04, 0.235])  # [h_placement, v_placement, h_size, v_size]
#cb = f.colorbar(ax2, cax=cbar_ax2, format='%.1f', extend='both') # double-edged colorbar
#cb.set_label('%s' % (cb_label), fontsize=16)                    # colorbar label

#f.subplots_adjust(right=0.8)
#cbar_ax3 = f.add_axes([0.85, 0.095, 0.04, 0.235])  # [h_placement, v_placement, h_size, v_size]
#cb = f.colorbar(ax3, cax=cbar_ax3, format='%.1f', extend='both') # double-edged colorbar
#cb.set_label('%s' % (cb_label), fontsize=16)                    # colorbar label

plt.savefig("%sCDOD_latvsLs_dsrunvsrefrun.png" % (fpath), bbox_inches='tight') 
exit()
## Temperature PLOT
temp_t = np.matrix.transpose(temp_avg_t)
temp_t = temp_t[:,l_1:l_2]

fig = plt.figure(figsize=(10,10), dpi=100)
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(t_t,lat_t,temp_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
plt.xlabel('Solar longitude / degrees', fontsize=14, labelpad=10)
plt.ylabel('Latitude / degrees', fontsize=14, labelpad=10)

# Extra Markers
ax.plot(c[0,0],c[0,1],'o',color='y',markersize=10)
 
# Ticks
ax.set_xticks(major_ticksx)                                                       
ax.set_xticks(minor_ticksx, minor=True)                                       
ax.set_yticks(major_ticksy)                                                       
ax.set_yticks(minor_ticksy, minor=True)
ax.tick_params(axis='both', labelsize=12, pad=10)

# Colour bar
cb = plt.colorbar(format='%.2f', extend='both')
cb.set_label('Temperature difference / K')
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator

plt.axis('tight')
plt.savefig("%sSurfTempDiff_LatvsTime_FY_uavg_tavg.png" % (fpath), bbox_inches='tight')

## Surface pressure and Atmospheric density at surface PLOT
ps_t = np.matrix.transpose(pres_avg_t)
ps_t = ps_t[:,l_1:l_2]
rho_t = np.matrix.transpose(rho_avg_t)
rho_t = rho_t[:,l_1:l_2]

f, axarr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t_t
y = lat_t
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'

cb_label = 'Atmospheric pressure difference / Pa'
cb_label2 = 'Atmospheric density difference / kg / $m^3$'
 
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')

ax1 = axarr[0].pcolormesh(x, y, ps_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
axarr[0].axis('tight') 
axarr[0].plot(c[0,0],c[0,1],'o',color='y',markersize=10)
axarr[0].set_xticks(major_ticksx)
axarr[0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0].set_yticks(major_ticksy)                                                       
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('(a)', fontsize=18)
axarr[0].tick_params(axis='both', labelsize=14)

ax2 = axarr[1].pcolormesh(x, y, rho_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
axarr[1].plot(c[0,0],c[0,1],'o',color='y',markersize=10)
axarr[1].set_title('(b)', fontsize=18)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.54, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax1, cax=cbar_ax, format='%.0f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=14)                    # colorbar label

cbar_ax2 = f.add_axes([0.85, 0.1, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb2 = f.colorbar(ax2, cax=cbar_ax2, format='%.1e', extend='both') # double-edged colorbar
cb2.set_label('%s' % (cb_label2), fontsize=14)                    # colorbar label

plt.savefig("%sPresDensDiff_LatvsLs_zonavg_tavg.png" % (fpath), bbox_inches='tight') 

# Zonal wind PLOT
u_t = np.matrix.transpose(u_avg_t)
u_t = u_t[:,l_1:l_2]
fig = plt.figure( figsize=(10,10), dpi=100)
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(t_t,lat_t,u_t,norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
plt.xlabel('Solar longitude / degrees',fontsize=16)
plt.ylabel('Latitude / degrees',fontsize=16)
ax.plot(c[0,0],c[0,1],'o',color='y',markersize=10) 

ax.set_xticks(major_ticksx)
ax.set_xticks(minor_ticksx, minor=True)
ax.set_yticks(major_ticksy)
ax.set_yticks(minor_ticksy, minor=True)
 
cb = plt.colorbar(format='%.1f', extend='both')
cb.set_label('Zonal wind velocity difference / m / s')
tick_locator = ticker.MaxNLocator(nbins=7)
cb.locator = tick_locator
cb.update_ticks()

plt.axis('tight')
plt.savefig("%sZonalWindDiff_LatvsTime_FY_uavg_tavg.png" % (fpath), bbox_inches='tight')

# ALL FLUXES on one plot
fslw_t = np.matrix.transpose(fslw_avg_t[l_1:l_2,:]) # Incoming (surf) long wave (IR) radiation
ftlw_t = np.matrix.transpose(ftlw_avg_t[l_1:l_2,:]) # Outgoing (top) long wave (IR) radiation
fssw_t = np.matrix.transpose(fssw_avg_t[l_1:l_2,:]) # Incoming (surf) short wave (VL) radiation
ftsw_t = np.matrix.transpose(ftsw_avg_t[l_1:l_2,:]) # Outgoing (top) short wave (VL) radiation

f, axarr = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t_t
y = lat_t
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'
cb_label = 'Radiative flux difference / W / $m^2$' 
    
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')
 
ax1 = axarr[0,0].pcolormesh(x, y, fslw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r', vmax=10,vmin=-5)# vmin=np.min((np.min(fslw_t),np.min(ftlw_t),np.min(fssw_t),np.min(ftsw_t))), vmax=np.max((np.max(fslw_t),np.max(ftlw_t),np.max(fssw_t),np.max(ftsw_t))))
axarr[0,0].axis('tight')
axarr[0,0].plot(c[0,0],c[0,1],'o',color='y',markersize=10)
axarr[0,0].set_xticks(major_ticksx)
axarr[0,0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0,0].set_yticks(major_ticksy)                                                       
axarr[0,0].set_yticks(minor_ticksy, minor=True)
axarr[0,0].set_title('Incident flux at surface (LW) (a)', fontsize=10)
axarr[0,0].tick_params(axis='both', labelsize=10)

dv1 = make_axes_locatable(axarr[0,0])
cax1 = dv1.append_axes("right",size="5%",pad=0.05)
cb = f.colorbar(ax1,cax=cax1, format='%.1f', extend='both')
cb.set_label('%s' % (cb_label), fontsize=10)

ax2 = axarr[0,1].pcolormesh(x, y, ftlw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
axarr[0,1].plot(c[0,0],c[0,1],'o',color='y',markersize=10)
axarr[0,1].set_title('Outgoing flux at top (LW) (b)', fontsize=10)
axarr[0,1].tick_params(axis='both', labelsize=10)

dv2 = make_axes_locatable(axarr[0,1])
cax2 = dv2.append_axes("right",size="5%",pad=0.05)
cb2 = f.colorbar(ax2,cax=cax2, format='%.1f', extend='both')
cb2.set_label('%s' % (cb_label), fontsize=10)

ax3 = axarr[1,0].pcolormesh(x, y, fssw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r',vmax=5,vmin=-30)
axarr[1,0].plot(c[0,0],c[0,1],'o',color='y',markersize=10)
axarr[1,0].set_title('Incident flux at surface (SW) (c)', fontsize=10)
axarr[1,0].tick_params(axis='both', labelsize=10)

dv3 = make_axes_locatable(axarr[1,0])
cax3 = dv3.append_axes("right",size="5%",pad=0.05)
cb3 = f.colorbar(ax3,cax=cax3, format='%.1f', extend='both')
cb3.set_label('%s' % (cb_label), fontsize=10)

ax4 = axarr[1,1].pcolormesh(x, y, ftsw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
axarr[1,1].plot(c[0,0],c[0,1],'o',color='y',markersize=10)
axarr[1,1].set_title('Outgoing flux at top (SW) (d)', fontsize=10)
axarr[1,1].tick_params(axis='both', labelsize=10)

dv4 = make_axes_locatable(axarr[1,1])
cax4 = dv4.append_axes("right",size="5%",pad=0.05)
cb4 = f.colorbar(ax4,cax=cax4, format='%.1f', extend='both')
cb4.set_label('%s' % (cb_label), fontsize=10)   

# Colorbar creation and placement
#f.subplots_adjust(right=0.8)
#cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8])  # [h_placement, v_placement, h_size, v_size]
#cb = f.colorbar(ax3, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
#cb.set_label('%s' % (cb_label), fontsize=14)                    # colorbar label

plt.savefig("%sfluxes_latvsLs_zonavg_tavg.png" % (fpath), bbox_inches='tight')

### Short-term Temperature and Heating rates at exact location vs Altitude (put in particle size or mmr)

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

### Plot explaination
# Storm starts at tstep=96, which is midnight of sol 8 relative to (0,0). However the storm is at 135W (midpoint). 
# So 360/24 = 15deg for each hour, meaning local time at 135W is 135/15=9hrs behind (0,0) local time, so at dust storm time insertion it is 15:00 locally. 
# We want the plot to start 1 day before the storm, which will be at tstep=84, since each tstep accounts for 2 hours.
# From tstep=84 we push forward 10 hours for approximate midnight plot start
### In reality the plot starts at 01:00 the night before the storm, local time 135W.
f, axarr = plt.subplots(3, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)

tl1, tl2 =89, 125
al=15
latl, lonl=6, 9

d1 = tempa[1][tl1:tl2,:al,latl,lonl]
d2 = tempb[1][tl1:tl2,:al,latl,lonl]
d3 = (88800/24.)*(sw_htrta[1][tl1:tl2,:al,latl,lonl] + lw_htrta[1][tl1:tl2,:al,latl,lonl]) # NET heat rate (SW cooling, LW heating) changed from K/s to K/hr
d4 = rdusta[1][tl1:tl2,:al,latl,lonl]
d5 = dustqa[1][tl1:tl2,:al,latl,lonl]

data  = np.matrix.transpose(d2)
data2 = np.matrix.transpose(d1)
data3 = np.matrix.transpose(d3) 
data4 = np.matrix.transpose(d4)
data5 = np.matrix.transpose(d5)

y = alt[:al]
x = np.linspace(0,72,data.shape[1])
xlabel = 'Local time / hr'
ylabel = 'Altitude above surface / km'
cb_label = 'Temperature / K' 
cb_label2 = 'Net heating rate / K / hr'

major_ticksx = np.arange(0,np.max(x)+1,6)
minor_ticksx = np.arange(0,np.max(x),2)
major_ticksy = np.arange(0,np.max(y)+1,10)
minor_ticksy = np.arange(0,np.max(y),2)

# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=16, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=16, va='center', rotation='vertical')

ax1 = axarr[0].pcolormesh(x, y, data, cmap='jet')
axarr[0].axis('tight')
axarr[0].set_xticks(major_ticksx)
axarr[0].set_yticks(major_ticksy)
axarr[0].set_xticks(minor_ticksx, minor=True)
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('Reference run (a)', fontsize=10)
axarr[0].tick_params(axis='both', labelsize=14)

ax2 = axarr[1].pcolormesh(x, y, data2, cmap='jet')
axarr[1].set_title('Dust storm run (b)', fontsize=10)
axarr[1].tick_params(axis='both', labelsize=14)
axarr[1].add_patch(mpl.patches.Rectangle((14, 0), 24, 9, facecolor="none", linestyle='dashed'))

ax3 = axarr[2].pcolormesh(x, y, data3, cmap='RdBu_r', vmax=10, vmin=-10)
axarr[2].set_title('Dust storm run (c)', fontsize=10)
axarr[2].tick_params(axis='both', labelsize=14)
axarr[2].add_patch(mpl.patches.Rectangle((14, 0), 24, 9, facecolor="none", linestyle='dashed'))

lvl = np.array([10**-6,10**-5,1*10**-4,10**-3]) # Contour levels
ax = axarr[2].contour(x,y,data5,lvl,colors='k',linewidth=3,locator=ticker.LogLocator())
plt.clabel(ax, fontsize=9, inline=1,fmt='%2.0e')

f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.38, 0.02, 0.52])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax1, cax=cbar_ax, format='%.0f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=16)                    # colorbar label

f.subplots_adjust(right=0.8)
cbar_ax2 = f.add_axes([0.85, 0.095, 0.02, 0.235])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax3, cax=cbar_ax2, format='%.0f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label2), fontsize=16)                    # colorbar label

#locs,labels = py.xticks()
#py.xticks(locs, map(lambda x: "%02d" % x, locs*1e9))

plt.savefig("%sheating.png" % (fpath), bbox_inches='tight')
plt.close('all')

### Time series plots
# settings
s_l = [-2.05, -6.12, 242.7]  # landing site marking on plot (actually for 244.7, Ls is messed up)
ticky_latlon = [60,10,30,10]      # tick settings [xmajor,xminor,ymajor,yminor] ticks
ticky_latalt = [60,10,20,10]
int_Ls = int(np.ceil(Ls.shape[0]/(9*Months)))

# Dust particle size contours
rd_ds1 = {}
rd_ds1[0] = alt[:alt_half]
rd_ds1[1] = lat_t
rd_ds1[2] = rd_avg[:,:alt_half,:]

# dust mmr average difference contours
dqd_ds = {}
dqd_ds[0] = alt[:alt_half]
dqd_ds[1] = lat_t
dqd_ds[2] = dqd_avg[:,:alt_half,:]

wind = {}
wind[0] = u_diff
wind[1] = v_diff 

## Dust storm 1 Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sDustqdiff_latlon_tseries_ds1.png' % (fpath), mola)

## Dust storm time series temp (time,lat,lon)
plt_timeseries(temp_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Temperature difference / K', int_Ls, '%stempdiff_latlon_tseries_ds1.png' % (fpath), mola)
## Dust storm time series surftemp (time,laplt_timeseries(dustq_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sDustqdiff_latlon_tseries_ds1.png' % (fpath), mola)t,lon)
plt_timeseries(tsurf_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Surface temperature difference / K', int_Ls, '%ssurftempdiff_latlon_tseries_ds1.png' % (fpath), mola)
## Dust storm time series pressure (time,lat,lon)
plt_timeseries(rho_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Atmospheric density difference / kg / $m^3$', int_Ls, '%sdensdiff_latlon_tseries_ds1.png' % (fpath), mola)
## 
plt_timeseries(pres_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Atmospheric pressure difference / Pa', int_Ls, '%spresdiff_latlon_tseries_ds1.png' % (fpath), mola)
## 
plt_timeseries(ps_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Surface pressure difference / Pa', int_Ls, '%ssurfpresdiff_latlon_tseries_ds1.png' % (fpath), mola)
## 
plt_timeseries(u_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Zonal wind velocity / m / s', int_Ls, '%sudiff_latlon_tseries_ds1.png' % (fpath), mola)
## fluxes
ftnet_diff = ftsw_diff + ftlw_diff
fsnet_diff = fssw_diff + fslw_diff

plt_timeseries(ftnet_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Net flux at top difference / W / $m^2$', int_Ls, '%snettop_latlon_tseries_ds1.png' % (fpath), mola)
plt_timeseries(fsnet_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Net flux at surface difference / W / $m^2$', int_Ls, '%snetsurf_latlon_tseries_ds1.png' % (fpath), mola)

plt_timeseries(ftlw_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Outgoing flux at top (LW) difference / W / $m^2$', int_Ls, '%sftlw_latlon_tseries_ds1.png' % (fpath), mola)
plt_timeseries(ftsw_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Outgoing flux at top (SW) difference / W / $m^2$', int_Ls, '%sftsw_latlon_tseries_ds1.png' % (fpath), mola)
plt_timeseries(fslw_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Incident flux at surface (LW) difference / W / $m^2$', int_Ls, '%sfslw_latlon_tseries_ds1.png' % (fpath), mola)
plt_timeseries(fssw_diff[l_1:l_2,:,:], lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Incident flux at surface (SW) difference / W / $m^2$', int_Ls, '%sfssw_latlon_tseries_ds1.png' % (fpath), mola)

## Dust lifting and sedimentation rates (changed from kg / $m^2$ s to kg / $m^2$ hr)
plt_timeseries(dqssed_diff[l_1:l_2,:,:]*(88800/24.), lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust sedimentation difference / kg / $m^2$ hr', int_Ls, '%sdqssed_latlon_tseries_ds1.png' % (fpath), mola)
plt_timeseries(dqsdev_diff[l_1:l_2,:,:]*(88800/24.), lon_t, lat_t, Ls, 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust devil lifting amount difference / kg / $m^2$ hr', int_Ls, '%sdqsdev_latlon_tseries_ds1.png' % (fpath), mola)

alt_t = alt # Height of 20.9km
dustq_diff_altlon = dustqa[1][l_1:l_2,:,8,:] - dustqb[1][l_1:l_2,:,8,:]
temp_diff_altlon = tempa[1][l_1:l_2,:,8,:] - tempb[1][l_1:l_2,:,8,:]

plt_timeseries(temp_diff_altlon, lon_t, alt_t, Ls, 4,4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Temperature difference / K', int_Ls, '%stemp_altlon_tseries_ds1.png' % (fpath))

plt_timeseries(dustq_diff_altlon, lon_t, alt_t, Ls, 4,4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sdustq_altlon_tseries_ds1.png' % (fpath))

plt.close('all')
#fix trigger to change from contour to pcolormesh
#fix trigger to label plots properly
#fix trigger to change the magnitude of units in colorbar
