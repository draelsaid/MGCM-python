# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016

import matplotlib as mpl
#mpl.use('Agg') # removes need for X-Server (graphics in linux). For qsub only.

import numpy as np
import matplotlib.colors as colors
import matplotlib.pyplot as plt
import datetime

from mpl_toolkits.axes_grid1 import make_axes_locatable
from mars_time import MarsTime
from scipy.io import *
from matplotlib import cm,ticker
from matplotlib.ticker import FormatStrFormatter

# Homemade
from MMM_plt_tseries import *
from MidPointNorm import *
#from plt_quiver import *

# Prints EVERYTHING inside a variable without holding back (intended for diagnostic)
np.set_printoptions(threshold=np.inf)

# Abbreviate sol_ls conversion function
sol_Ls=MarsTime().sol_ls

# Moving average 
def moving_average(a, n=3) :
    ret = np.cumsum(a, dtype=float)
    ret[n:] = ret[n:] - ret[:-n]
    return ret[n - 1:] / n

# custom rounding
def myround(x, base=5):
    return int(base * round(float(x)/base))

# Grab topography from surface.nc or mola32.nc file
ml = netcdf.netcdf_file('/padata/mars/users/aes442/mgcm_data/surface.nc','r')

mola = {}

mola[0] = ml.variables['latitude'][:]
mola[1] = ml.variables['longitude'][:]
mola[2] = ml.variables['zMOL'][:]

# Initialise dictionaries - due to data size
Ls_m = {}
hgt, phtot = {}, {}

tempa, tempb = {}, {}
ptota, ptotb = {}, {}
ua, ub = {}, {}
va, vb = {}, {}
dusta, dustb = {}, {}

swdwa, swdwb = {}, {}
lwdwa, lwdwb = {}, {}
swupa, swupb = {}, {}
lwupa, lwupb = {}, {}
hfxuwa, hfxuwb={}, {}
taua, taub = {}, {}

# Number of files in run

# DATA SHAPE >>> 
# bottom_top = ALTITUDE
# south_north= LATITUDE
# west_east  = LONGITUDE
#                   4D: (Time, bottom_top, south_north, west_east)
#                   3D: (Time, south_north, west_east)

m=0
for i in xrange(2,3):
 m=m+1
 
 mgcm = "LMD_MMGCM"
 rundira = "3_ds_callphys"
 rundirb = "3_ref_callphys"
 filename = "wrfout_d01_2024-01-0%i_02:00:00" % (i)
 
 a = netcdf.netcdf_file("/padata/mars/users/aes442/RUNS/R-%s/%s/%s" % (mgcm,rundira,filename),'r')
 b = netcdf.netcdf_file("/padata/mars/users/aes442/RUNS/R-%s/%s/%s" % (mgcm,rundirb,filename),'r')

## Referencing spatial/temporal variables
# Altitude
 znu      = a.variables['ZNU'][:]        
 znw      = a.variables['ZNW'][:]
 hgt[m]   = a.variables['HGT'][:]
 phtot[m] = a.variables['PHTOT'][:]

# Latitude and Longitude
 xlat = a.variables['XLAT'][:][0,:,0]
 xlon = a.variables['XLONG'][0,0,:]

# Required variables 
 dusta[m]  = a.variables['QDUST'][:]
 tempa[m]  = a.variables['TEMP'][:]
 ptota[m] = a.variables['PTOT'][:]
 ua[m]     = a.variables['U'][:]
 va[m]     = a.variables['V'][:]
  
 swdwa[m]  = a.variables['SWDOWNZ'][:]
 lwdwa[m]  = a.variables['LWDOWNZ'][:]
 swupa[m]  = a.variables['SWUP'][:]
 lwupa[m]  = a.variables['LWUP'][:]
 hfxuwa[m] = a.variables['HFX'][:]
 taua[m] = a.variables['TAU'][:]
 
## b
 dustb[m]  = b.variables['QDUST'][:]
 tempb[m]  = b.variables['TEMP'][:]
 ptotb[m] = b.variables['PTOT'][:]
 ub[m]     = b.variables['U'][:]
 vb[m]     = b.variables['V'][:]
  
 swdwb[m]  = a.variables['SWDOWNZ'][:]
 lwdwb[m]  = a.variables['LWDOWNZ'][:]
 swupb[m]  = a.variables['SWUP'][:]
 lwupb[m]  = a.variables['LWUP'][:]
 hfxuwb[m] = a.variables['HFX'][:]
 taub[m] = b.variables['TAU'][:]

# Calculate approximate HEIGHT from geopotential
alt = np.zeros((phtot[1].shape[1]))
alt = phtot[1][0,:,0,0] / 3.72 - hgt[1][0,0,0]

print "Latitude: %i, (%.2f - %.2f) || Longitude: %i, (%.2f - %.2f) || Model levels: %i => Alt Min:%.3f | Alt Max:%.3f " % (xlat.shape[0],xlat[0],xlat[-1],xlon.shape[0],xlon[0],xlon[-1],alt.shape[0],alt[0],alt[-1])

# Time

fpath = "/home/physastro/aes442/results/MMM_dust/"

#########################################################################################

## DATA

## MOLA TOPOGRAPHY
dd1 = np.where(np.absolute(mola[0]-xlat[0])==np.min(np.absolute(mola[0]-xlat[0])))[0][0]
dd2 = np.where(np.absolute(mola[0]-xlat[-1])==np.min(np.absolute(mola[0]-xlat[-1])))[0][0]
dd3 = np.where(np.absolute(mola[1]-xlon[0])==np.min(np.absolute(mola[1]-xlon[0])))[0][0]
dd4 = np.where(np.absolute(mola[1]-xlon[-1])==np.min(np.absolute(mola[1]-xlon[-1])))[0][0]

topg = {}
topg[0] = mola[0][dd2-1:dd1+1] # lat
topg[1] = mola[1][dd3-1:dd4+1] # lon
topg[2] = mola[2][dd2-1:dd1+1,dd3-1:dd4+1] # (lat,lon) mola map

day = 1
hr = 13
alt = 0

## variable[day][hour, elevation, lat, lon]
## day average

#ut = ua[day][:,alt,0:100,0:100].sum(axis=0)/ua[day].shape[0] - ub[day][:,alt,0:100,0:100].sum(axis=0)/ub[day].shape[0]
#vt = va[day][:,alt,0:100,0:100].sum(axis=0)/va[day].shape[0] - vb[day][:,alt,0:100,0:100].sum(axis=0)/vb[day].shape[0]

#data = tempa[day][:,alt,:,:].sum(axis=0)/tempa[day].shape[0] - tempb[day][:,alt,:,:].sum(axis=0)/tempb[day].shape[0]
#data2= ptota[day][:,alt,:,:].sum(axis=0)/ptota[day].shape[0] - ptotb[day][:,alt,:,:].sum(axis=0)/ptotb[day].shape[0]

ut = ua[day][hr,alt,0:100,0:100] - ub[day][hr,alt,0:100,0:100]
vt = va[day][hr,alt,0:100,0:100] - vb[day][hr,alt,0:100,0:100]

data = tempa[day][hr,alt,:,:] - tempb[day][hr,alt,:,:]
data2= ptota[day][hr,alt,:,:] - ptotb[day][hr,alt,:,:]

## PLOTS

# Common settings (ticks)
t = np.arange(0,24)

# Longitude
major_ticksx = np.arange(np.floor(xlon[0]), np.ceil(xlon[-1]), 10)
minor_ticksx = np.arange(np.floor(xlon[0]), np.ceil(xlon[-1]), 5)

# Latitude
major_ticksy = np.arange(np.floor(xlat[0]), np.ceil(xlat[-1]), 10)
minor_ticksy = np.arange(np.floor(xlat[0]), np.ceil(xlat[-1]), 5)

## PLOT temperature/winds/topography
f, axarr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = xlon
y = xlat

xlabel = 'Longitude / degrees'
ylabel = 'Latitude / degrees'
cblabel= 'Temperature difference / K'
plt.xlabel(xlabel, fontsize=14, labelpad=10)
plt.ylabel(ylabel, fontsize=14, labelpad=10)

# Main plot
ax = axarr.pcolormesh(x, y, data, cmap='RdBu_r', norm=MidPointNorm(midpoint=0.))

# Secondary plot
ax2 = axarr.quiver(xlon, xlat, ut, vt, scale=2**2, units='y', width=0.1)
aq = axarr.quiverkey(ax2, 0.815, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')

# Topography
lvls = [-5,0,5,10,15]
ax3 = axarr.contour(topg[1], topg[0], topg[2], lvls, colors='k')
 
# Ticks
axarr.set_xticks(major_ticksx)                                                       
axarr.set_xticks(minor_ticksx, minor=True)                                       
axarr.set_yticks(major_ticksy)                                                       
axarr.set_yticks(minor_ticksy, minor=True)
axarr.tick_params(axis='both', labelsize=12, pad=10)

axarr.axis('tight')

# Colour bar
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8])    # [h_place, v_place, h_size, v_size]
cb = f.colorbar(ax, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cblabel), fontsize=16)                    # colorbar label

plt.axis('tight')
plt.savefig("%sMMM_temp_uvwind_mola_latvslon.png" % (fpath), bbox_inches='tight')
#plt.close('all')
exit()

## PLOT pressure/topography
f, axarr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = xlon
y = xlat

cblabel= 'Pressure difference / Pa'
plt.xlabel(xlabel, fontsize=14, labelpad=10)
plt.ylabel(ylabel, fontsize=14, labelpad=10)

# Main plot
ax = axarr.pcolormesh(x, y, data2, cmap='RdBu_r', norm=MidPointNorm(midpoint=0.))

# Topography
lvls = [-5,0,5,10,15]
ax3 = axarr.contour(topg[1], topg[0], topg[2], lvls, colors='k')
 
# Ticks
axarr.set_xticks(major_ticksx)                                                       
axarr.set_xticks(minor_ticksx, minor=True)                                       
axarr.set_yticks(major_ticksy)                                                       
axarr.set_yticks(minor_ticksy, minor=True)
axarr.tick_params(axis='both', labelsize=12, pad=10)

axarr.axis('tight')

# Colour bar
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8])    # [h_place, v_place, h_size, v_size]
cb = f.colorbar(ax, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cblabel), fontsize=16)                    # colorbar label

plt.axis('tight')
plt.savefig("%sMMM_pressure_latvslon.png" % (fpath), bbox_inches='tight')
plt.close('all')



# plot
f, axarr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = xlat
y = alt[:50]/1000.

# DATA temp alt/lon
aa = tempa[day][15,:,:,:]
bb = tempb[day][15,:,:,:]

# time averaging
#aa = tempa[day].sum(axis=0)/tempa[day].shape[0]
#bb = tempb[day].sum(axis=0)/tempb[day].shape[0]

# lon averaging
aa2 = aa.sum(axis=2)/aa.shape[2]
bb2 = bb.sum(axis=2)/bb.shape[2]

tmp_data = aa2 - bb2

# Labels
xlabel = 'Latitude / degrees'
ylabel = 'Altitude / km'
cb_label = 'Temperature difference / K' 
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')

# Altitude
major_ticksy = np.arange(np.floor(alt[0]), np.ceil(alt[-1]), 10)
minor_ticksy = np.arange(np.floor(alt[0]), np.ceil(alt[-1]), 5)

# Latitude
major_ticksx = np.arange(np.floor(xlat[0]), np.ceil(xlat[-1]), 10)
minor_ticksx = np.arange(np.floor(xlat[0]), np.ceil(xlat[-1]), 5)

ax = axarr.pcolormesh(x, y, tmp_data, cmap='RdBu_r', norm=MidPointNorm(midpoint=0.))

axarr.axis('tight')

# Colour bar
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8])    # [h_place, v_place, h_size, v_size]
cb = f.colorbar(ax, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=16)                    # colorbar label

#plt.axis('tight')
plt.savefig("%sMMM_tmp_altvslat.png" % (fpath), bbox_inches='tight')
plt.close('all')


## PLOT radiative fluxes (zonally averaged)

# data
data1 = hfxuwa[day].sum(axis=2)/hfxuwa[day].shape[2] - hfxuwb[day].sum(axis=2)/hfxuwb[day].shape[2]
data2 = lwdwa[day].sum(axis=2)/lwdwa[day].shape[2] - lwdwb[day].sum(axis=2)/lwdwb[day].shape[2]
data3 = swdwa[day].sum(axis=2)/swdwa[day].shape[2] - swdwb[day].sum(axis=2)/swdwb[day].shape[2]
data1 = data1.T
data2 = data2.T
data3 = data3.T

# plot 
f, axarr = plt.subplots(1, 2, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t
y = xlat

# Labels
xlabel = 'Time / hours'
ylabel = 'Latitude / degrees'
cb_label = 'Radiative flux difference / W / $m^2$' 
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')

# Time (hours)
major_ticksx = np.arange(0, 24, 2)
minor_ticksx = np.arange(0, 24, 0.5)

# Latitude
major_ticksy = np.arange(np.floor(xlat[0]), np.ceil(xlat[-1]), 10)
minor_ticksy = np.arange(np.floor(xlat[0]), np.ceil(xlat[-1]), 5)

ax1 = axarr[0].pcolormesh(x, y, data2, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')# 
axarr[0].axis('tight')
axarr[0].set_xticks(major_ticksx)
axarr[0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0].set_yticks(major_ticksy)                                                       
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('Outgoing heat flux at surface (a)', fontsize=10)
axarr[0].tick_params(axis='both', labelsize=10)

dv1 = make_axes_locatable(axarr[0])
cax1 = dv1.append_axes("right",size="5%",pad=0.05)
cb = f.colorbar(ax1,cax=cax1, format='%.1f', extend='both')
cb.set_label('%s' % (cb_label), fontsize=10)

ax2 = axarr[1].pcolormesh(x, y, data3, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
axarr[1].set_title('Incident flux at surface (SW) (b)', fontsize=10)
axarr[1].tick_params(axis='both', labelsize=10)

dv2 = make_axes_locatable(axarr[1])
cax2 = dv2.append_axes("right",size="5%",pad=0.05)
cb2 = f.colorbar(ax2,cax=cax2, format='%.1f', extend='both')
cb2.set_label('%s' % (cb_label), fontsize=10)

plt.axis('tight')
plt.savefig("%sMMM_influxes_latvstime.png" % (fpath), bbox_inches='tight')
plt.close('all')

# SENSIBLE heat flux
f, axarr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t
y = xlat

ax1 = axarr.pcolormesh(x, y, data1, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
axarr.axis('tight')
axarr.set_xticks(major_ticksx)
axarr.set_xticks(minor_ticksx, minor=True)                                           
axarr.set_yticks(major_ticksy)                                                       
axarr.set_yticks(minor_ticksy, minor=True)
axarr.set_title('Outgoing sensible heat flux at surface', fontsize=10)
axarr.tick_params(axis='both', labelsize=10)

dv1 = make_axes_locatable(axarr)
cax1 = dv1.append_axes("right",size="5%",pad=0.05)
cb = f.colorbar(ax1,cax=cax1, format='%.1f', extend='both')
cb.set_label('%s' % (cb_label), fontsize=10)

plt.axis('tight')
plt.savefig("%sMMM_outflux_latvstime.png" % (fpath), bbox_inches='tight')
plt.close('all')

## Dust storm 1 Time series dustq (mmr) (time, lat, lon)

d_data = dusta[day][:,0,0:100,0:100] - dustb[day][:,0,0:100,0:100]
tau_data = dusta[day][:,0,0:100,0:100] - dustb[day][:,0,0:100,0:100]

plt_timeseries(d_data, xlon, xlat, t, 4,6, 'Longitude / degrees', 'Latitude / degrees', 'Hour: ', 'Dust MMR difference / kg / kg', 1, '%sMMM_dustqdiff_latlon_tseries.png' % (fpath), topg)

plt_timeseries(d_data, xlon, xlat, t, 4,6, 'Longitude / degrees', 'Latitude / degrees', 'Hour: ', 'Dust MMR difference / kg / kg', 1, '%sMMM_dustqdiff_latlon_tseries.png' % (fpath), topg)

plt.close('all')

