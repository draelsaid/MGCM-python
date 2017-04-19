# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016

import matplotlib as mpl
#mpl.use('Agg') # removes need for X-Server (graphics in linux). For qsub only.

import numpy as np
import pylab as py
import matplotlib.colors as colors
import matplotlib.pyplot as plt

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

# Grab topography from surface.nc or mola32.nc file
ml = netcdf.netcdf_file('/padata/mars/users/aes442/mgcm_data/surface.nc','r')

mola = {}

mola[0] = ml.variables['latitude'][:]
mola[1] = ml.variables['longitude'][:]
mola[2] = ml.variables['zMOL'][:]

# Import data from Luca's TES dust files for comparison
a = netcdf.netcdf_file('/padata/mars/users/aes442/mgcm_data/dust_MY28.nc','r')

d_lat_s = a.variables['latitude'][:]
d_lon_s = a.variables['longitude'][:]
d_t   = a.variables['Time'][:]
d_d   = a.variables['dustop'][:]

d_lat = np.linspace(-90,90,d_lat_s.shape[0])
d_lon = np.linspace(-180,180,d_lon_s.shape[0])

# Number of months in comparison (always add 1 because of Python indexing)
Months = 2   # No. of months
amth = 30     # Actual month 

# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,Months):
 mgcm = "MGCM_v5-1"
 rundira = "new_ds"
 rundirb = "new_ref"
 month = ("m%s" % (amth)) # CHANGE
 filename = "diagfi.nc"
 
 a = netcdf.netcdf_file("/padata/alpha/users/aes442/RUNS/R-%s/%s/%s/%s" % (mgcm,rundira,month,filename),'r')
 b = netcdf.netcdf_file("/padata/alpha/users/aes442/RUNS/R-%s/%s/%s/%s" % (mgcm,rundirb,month,filename),'r')
 
 lat     = a.variables['lat'][:]
 lon     = a.variables['lon'][:]
 sigma   = a.variables['sigma'][:]
 t_m     = a.variables['time'][:]
 Ls_m[i] = a.variables['Ls'][:]
 
 psa[i]         = a.variables['ps'][:]
 presa[i]       = a.variables['pressure'][:]
 tempa[i]       = a.variables['temp'][:]
 tsurfa[i]      = a.variables['tsurf'][:]
 ua[i]          = a.variables['u'][:]
 va[i]          = a.variables['v'][:]
 dustqa[i]      = a.variables['dustq'][:]
 dustNa[i]      = a.variables['dustN'][:]
 rhoa[i]        = a.variables['rho'][:]
 fluxsurflwa[i] = a.variables['fluxsurf_lw'][:]
 fluxsurfswa[i] = a.variables['fluxsurf_sw'][:]
 fluxtoplwa[i]  = a.variables['fluxtop_lw'][:]
 fluxtopswa[i]  = a.variables['fluxtop_sw'][:]
 taua[i]        = a.variables['taudustvis'][:]
 rdusta[i]      = a.variables['reffdust'][:]
 lw_htrta[i]    = a.variables['lw_htrt'][:]
 sw_htrta[i]    = a.variables['sw_htrt'][:]
 dqsseda[i]     = a.variables['dqssed'][:]
 dqsdeva[i]     = a.variables['dqsdev'][:]
 
 psb[i]         = b.variables['ps'][:]
 presb[i]       = b.variables['pressure'][:]
 tempb[i]       = b.variables['temp'][:]
 tsurfb[i]      = b.variables['tsurf'][:]
 ub[i]          = b.variables['u'][:]
 vb[i]          = b.variables['v'][:]
 dustqb[i]      = b.variables['dustq'][:]
 dustNb[i]      = b.variables['dustN'][:]
 rhob[i]        = b.variables['rho'][:]
 fluxsurflwb[i] = b.variables['fluxsurf_lw'][:]
 fluxsurfswb[i] = b.variables['fluxsurf_sw'][:]
 fluxtoplwb[i]  = b.variables['fluxtop_lw'][:]
 fluxtopswb[i]  = b.variables['fluxtop_sw'][:]
 taub[i]        = b.variables['taudustvis'][:]
 rdustb[i]      = b.variables['reffdust'][:]
 lw_htrtb[i]    = b.variables['lw_htrt'][:]
 sw_htrtb[i]    = b.variables['sw_htrt'][:]
 dqssedb[i]     = b.variables['dqssed'][:]
 dqsdevb[i]     = b.variables['dqsdev'][:] 

# Calculate approximate HEIGHT from sigma (km)
alt = np.zeros((sigma.shape[0]))
for i in xrange(len(sigma)):
 alt[i] = -10.8*np.log(sigma[i])

print "Latitude: %i || Longitude: %i || Model levels: %i => Alt Min:%.3f | Alt Max:%.3f | Alt half: %.3f " % (lat.shape[0],lon.shape[0],sigma.shape[0],alt[0],alt[-1],alt[18])
alt_half=18 # 47.8km

# Get time dimension length
n = 0
for i in xrange(1,len(psa)+1,1): # len(psa) gives the number of months
 n = n + len(dustqa[i])          # len(dustqa[i]) gives the number of time steps in each month.
print ("Total time steps: %i" % (n))

## Ls vector
Ls_s = (Months-1)*30 # Number of solar longitudes for time vector for comparison
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

## Create all other variables, with altitude dimension removed
ps_a, ps_b           = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
temp_a, temp_b       = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
tsurf_a, tsurf_b     = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
u_a, u_b             = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
v_a, v_b             = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustq_a, dustq_b     = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustN_a, dustN_b     = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
rho_a, rho_b         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fslwa, fslwb         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fsswa, fsswb         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
ftlwa, ftlwb         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
ftswa, ftswb         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
tau_a, tau_b         = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
rdust_a, rdust_b     = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
lw_htrt_a, lw_htrt_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
sw_htrt_a, sw_htrt_b = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
pres_a, pres_b       = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dqssed_a, dqssed_b   = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dqsdev_a, dqsdev_b   = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))

# 3D Vars
ps_a, ps_b = psa[1][:,:,:], psb[1][:,:,:]
fslwa, fslwb = fluxsurflwa[1][:,:,:], fluxsurflwb[1][:,:,:]
fsswa, fsswb = fluxsurfswa[1][:,:,:], fluxsurfswb[1][:,:,:]
ftlwa, ftlwb = fluxtoplwa[1][:,:,:], fluxtoplwb[1][:,:,:]
ftswa, ftswb = fluxtopswa[1][:,:,:], fluxtopswb[1][:,:,:]
tau_a, tau_b = taua[1][:,:,:], taub[1][:,:,:]
tsurf_a, tsurf_b = tsurfa[1][:,:,:], tsurfb[1][:,:,:]
dqssed_a, dqssed_b = dqsseda[1][:,:,:], dqssedb[1][:,:,:]
dqsdev_a, dqsdev_b = dqsdeva[1][:,:,:], dqsdevb[1][:,:,:]

# 4D Vars
temp_a, temp_b   = tempa[1][:,1,:,:], tempb[1][:,1,:,:]
u_a, u_b   = ua[1][:,1,:,:], ub[1][:,1,:,:]
v_a, v_b   = va[1][:,1,:,:], vb[1][:,1,:,:]
dustq_a, dustq_b   = dustqa[1][:,1,:,:], dustqb[1][:,1,:,:]
dustN_a, dustN_b   = dustNa[1][:,1,:,:], dustNb[1][:,1,:,:]
rho_a, rho_b   = rhoa[1][:,1,:,:], rhob[1][:,1,:,:]
rdust_a, rdust_b = rdusta[1][:,1,:,:], rdustb[1][:,1,:,:]
lw_htrt_a, lw_htrt_b = lw_htrta[1][:,1,:,:], lw_htrtb[1][:,1,:,:]
sw_htrt_a, sw_htrt_b = sw_htrta[1][:,1,:,:], sw_htrtb[1][:,1,:,:]
pres_a, pres_b = presa[1][:,1,:,:], presb[1][:,1,:,:]

# Longitudal averaging

# Variables without longitude
temp_aa, temp_bb       = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
tsurf_aa, tsurf_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
u_aa, u_bb             = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustq_aa, dustq_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustN_aa, dustN_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rho_aa, rho_bb         = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rdust_aa, rdust_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
lw_htrt_aa, lw_htrt_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
sw_htrt_aa, sw_htrt_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
pres_aa, pres_bb       = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))

# 4D Vars
temp_aa, temp_bb = np.sum(tempa[1],axis=3)/tempa[1].shape[3], np.sum(tempb[1],axis=3)/tempb[1].shape[3]
u_aa, u_bb = np.sum(ua[1],axis=3)/ua[1].shape[3], np.sum(ub[1],axis=3)/ub[1].shape[3]
dustq_aa, dustq_bb = np.sum(dustqa[1],axis=3)/dustqa[1].shape[3], np.sum(dustqb[1],axis=3)/dustqb[1].shape[3]
dustN_aa, dustN_bb = np.sum(dustNa[1],axis=3)/dustNa[1].shape[3], np.sum(dustNb[1],axis=3)/dustNb[1].shape[3]
rho_aa, rho_bb = np.sum(rhoa[1],axis=3)/rhoa[1].shape[3], np.sum(rhob[1],axis=3)/rhob[1].shape[3]
rdust_aa, rdust_bb = np.sum(rdusta[1],axis=3)/rdusta[1].shape[3], np.sum(rdustb[1],axis=3)/rdustb[1].shape[3]
lw_htrt_aa, lw_htrt_bb = np.sum(lw_htrta[1],axis=3)/lw_htrta[1].shape[3], np.sum(lw_htrtb[1],axis=3)/lw_htrtb[1].shape[3]
sw_htrt_aa, sw_htrt_bb = np.sum(sw_htrta[1],axis=3)/sw_htrta[1].shape[3], np.sum(sw_htrtb[1],axis=3)/sw_htrtb[1].shape[3]
pres_aa, pres_bb = np.sum(presa[1],axis=3)/presa[1].shape[3], np.sum(presb[1],axis=3)/presb[1].shape[3]

# Calculate differences
dustq_diff = dustq_a - dustq_b
dustN_diff = dustN_a - dustN_b
temp_diff = temp_a - temp_b 
tsurf_diff = tsurf_a - tsurf_b 
ps_diff = ps_a - ps_b
rho_diff = rho_a - rho_b
u_diff = u_a - u_b
v_diff = v_a - v_b
rdust_diff = rdust_a - rdust_b
lw_htrt_diff = lw_htrt_a - lw_htrt_b
sw_htrt_diff = sw_htrt_a - sw_htrt_b
pres_diff = pres_a - pres_b
dqssed_diff = dqssed_a - dqssed_b
dqsdev_diff = dqsdev_a - dqsdev_b

fslw_diff = fslwa - fslwb
fssw_diff = fsswa - fsswb
ftlw_diff = ftlwa - ftlwb
ftsw_diff = ftswa - ftswb
 
t_d = temp_aa - temp_bb
pres_d = pres_aa - pres_bb
ts_d = tsurf_aa - tsurf_bb
dq_d = dustq_aa - dustq_bb
dN_d = dustN_aa - dustN_bb
rho_d = rho_aa - rho_bb
u_d = u_aa - u_bb
rdust_d = rdust_aa - rdust_bb
lw_htrt_d = lw_htrt_aa - lw_htrt_bb
sw_htrt_d = sw_htrt_aa - sw_htrt_bb

# Zonal averaging (time,lat)
temp_avg = np.sum(temp_a,axis=2)/temp_a.shape[2] - np.sum(temp_b,axis=2)/temp_b.shape[2]
tsurf_avg = np.sum(tsurf_a,axis=2)/tsurf_a.shape[2] - np.sum(tsurf_b,axis=2)/tsurf_b.shape[2]
ps_avg = np.sum(ps_a,axis=2)/ps_a.shape[2] - np.sum(ps_b,axis=2)/ps_b.shape[2]
pres_avg = np.sum(pres_a,axis=2)/pres_a.shape[2] - np.sum(pres_b,axis=2)/pres_b.shape[2]
u_avg = np.sum(u_a,axis=2)/u_a.shape[2] - np.sum(u_b,axis=2)/u_b.shape[2]
rho_avg = np.sum(rho_a,axis=2)/rho_a.shape[2] - np.sum(rho_b,axis=2)/rho_b.shape[2]
fssw_avg = np.sum(fsswa,axis=2)/fsswa.shape[2] - np.sum(fsswb,axis=2)/fsswb.shape[2]
fslw_avg = np.sum(fslwa,axis=2)/fslwa.shape[2] - np.sum(fslwb,axis=2)/fslwb.shape[2]
ftsw_avg = np.sum(ftswa,axis=2)/ftswa.shape[2] - np.sum(ftswb,axis=2)/ftswb.shape[2]
ftlw_avg = np.sum(ftlwa,axis=2)/ftlwa.shape[2] - np.sum(ftlwb,axis=2)/ftlwb.shape[2]
tau_a_avg = np.sum(tau_a,axis=2)/tau_a.shape[2]
tau_b_avg = np.sum(tau_b,axis=2)/tau_b.shape[2]
rdust_avg = np.sum(rdust_a,axis=2)/rdust_a.shape[2] - np.sum(rdust_b,axis=2)/rdust_b.shape[2]
lw_htrt_avg = np.sum(lw_htrt_a,axis=2)/lw_htrt_a.shape[2] - np.sum(lw_htrt_b,axis=2)/lw_htrt_b.shape[2]
sw_htrt_avg = np.sum(sw_htrt_a,axis=2)/sw_htrt_a.shape[2] - np.sum(sw_htrt_b,axis=2)/sw_htrt_b.shape[2]
 
temp_avg_ = np.sum(temp_b,axis=2)/temp_b.shape[2]
pres_avg_ = np.sum(pres_b,axis=2)/pres_b.shape[2]
tsurf_avg_ = np.sum(tsurf_b,axis=2)/tsurf_b.shape[2]
ps_avg_ = np.sum(ps_b,axis=2)/ps_b.shape[2]
u_avg_ = np.sum(u_b,axis=2)/u_b.shape[2]
rho_avg_ = np.sum(rho_b,axis=2)/rho_b.shape[2]
fssw_avg_ = np.sum(fsswb,axis=2)/fsswb.shape[2]
fslw_avg_ = np.sum(fslwb,axis=2)/fslwb.shape[2]
ftsw_avg_ = np.sum(ftswb,axis=2)/ftswb.shape[2]
ftlw_avg_ = np.sum(ftlwb,axis=2)/ftlwb.shape[2]

# from 35N to 55N Lat 
tmp_ = np.sum(np.sum(temp_avg_[:,7:11],axis=0)/n,axis=0)/4
tmps_ = np.sum(np.sum(tsurf_avg_[:,7:11],axis=0)/n,axis=0)/4
ps_ = np.sum(np.sum(ps_avg_[:,7:11],axis=0)/n,axis=0)/4
pres_ = np.sum(np.sum(pres_avg_[:,7:11],axis=0)/n,axis=0)/4
rho_ = np.sum(np.sum(rho_avg_[:,7:11],axis=0)/n,axis=0)/4
u_ = np.sum(np.sum(np.absolute(u_avg_[:,7:11]),axis=0)/n,axis=0)/4
fslw_ = np.sum(np.sum(fslw_avg_[:,7:11],axis=0)/n,axis=0)/4
fssw_ = np.sum(np.sum(fssw_avg_[:,7:11],axis=0)/n,axis=0)/4
ftlw_ = np.sum(np.sum(ftlw_avg_[:,7:11],axis=0)/n,axis=0)/4
ftsw_ = np.sum(np.sum(ftsw_avg_[:,7:11],axis=0)/n,axis=0)/4

tmp_1 = np.sum(np.sum(temp_avg[:,7:11],axis=0)/n,axis=0)/4
tmps_1 = np.sum(np.sum(tsurf_avg[:,7:11],axis=0)/n,axis=0)/4
ps_1 = np.sum(np.sum(ps_avg[:,7:11],axis=0)/n,axis=0)/4
pres_1 = np.sum(np.sum(pres_avg[:,7:11],axis=0)/n,axis=0)/4
rho_1 = np.sum(np.sum(rho_avg[:,7:11],axis=0)/n,axis=0)/4
u_1 = np.sum(np.sum(u_avg[:,7:11],axis=0)/n,axis=0)/4
fslw_1 = np.sum(np.sum(fslw_avg[:,7:11],axis=0)/n,axis=0)/4
fssw_1 = np.sum(np.sum(fssw_avg[:,7:11],axis=0)/n,axis=0)/4
ftlw_1 = np.sum(np.sum(ftlw_avg[:,7:11],axis=0)/n,axis=0)/4
ftsw_1 = np.sum(np.sum(ftsw_avg[:,7:11],axis=0)/n,axis=0)/4

print "AVERAGES:  tmp: %.2f || surf tmp: %.2f || press: %.2f || surf press: %.2f || dens: %.2f || zon wind: %.2f || fluxes (inLW: %.2f, outLW: %.2f, inSW: %.2f, outSW: %.2f). " % (tmp_, tmps_, pres_, ps_, rho_, u_, fslw_, ftlw_, fssw_, ftsw_)
#print tmp_1/tmp_, tmps_1/tmps_, pres_1/pres_, ps_1/ps_, rho_1/rho_, u_1/u_, fslw_1/fslw_, fssw_1/fssw_, ftlw_1/ftlw_, ftsw_1/ftsw_

# Time moving-point average of zonal average
nn=2 # Number of points to average over
t_avg = Ls[:-(nn-1)]

temp_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
pres_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
tsurf_avg_t   = np.zeros((t_avg.shape[0],lat.shape[0]))
ps_avg_t      = np.zeros((t_avg.shape[0],lat.shape[0]))
u_avg_t       = np.zeros((t_avg.shape[0],lat.shape[0]))
rho_avg_t     = np.zeros((t_avg.shape[0],lat.shape[0]))
fssw_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
fslw_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
ftsw_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
ftlw_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
rdust_avg_t   = np.zeros((t_avg.shape[0],lat.shape[0]))
lw_htrt_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))
sw_htrt_avg_t = np.zeros((t_avg.shape[0],lat.shape[0]))

for i in xrange(0,lat.shape[0]): 
 temp_avg_t[:,i]  = moving_average(temp_avg[:,i],n=nn)
 pres_avg_t[:,i]  = moving_average(pres_avg[:,i],n=nn)
 tsurf_avg_t[:,i] = moving_average(tsurf_avg[:,i],n=nn)
 ps_avg_t[:,i]    = moving_average(ps_avg[:,i],n=nn)
 u_avg_t[:,i]     = moving_average(u_avg[:,i],n=nn)
 rho_avg_t[:,i]   = moving_average(rho_avg[:,i],n=nn)
 fssw_avg_t[:,i]  = moving_average(fssw_avg[:,i],n=nn)
 fslw_avg_t[:,i]  = moving_average(fslw_avg[:,i],n=nn)
 ftsw_avg_t[:,i]  = moving_average(ftsw_avg[:,i],n=nn)
 ftlw_avg_t[:,i]  = moving_average(ftlw_avg[:,i],n=nn)
 rdust_avg_t[:,i] = moving_average(rdust_avg[:,i],n=nn)
 lw_htrt_avg_t[:,i] = moving_average(lw_htrt_avg[:,i],n=nn)
 sw_htrt_avg_t[:,i] = moving_average(sw_htrt_avg[:,i],n=nn)

############ TIME AVERAGE of differences ###################
nnn=nn
t_av = Ls[:-(nnn-1)]

td_avg    = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
pres_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
tds_avg   = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
dqd_avg   = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
dNd_avg   = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
rhod_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
ud_avg    = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
rd_avg    = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
lwhr_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
swhr_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))

for j in xrange(0,lat.shape[0],1):
 for i in xrange(0,sigma.shape[0],1):
  td_avg[:,i,j]   = moving_average(t_d[:,i,j],n=nnn)
  pres_avg[:,i,j] = moving_average(pres_d[:,i,j],n=nnn)
  tds_avg[:,i,j]  = moving_average(ts_d[:,i,j],n=nnn)
  dqd_avg[:,i,j]  = moving_average(dq_d[:,i,j],n=nnn)
  dNd_avg[:,i,j]  = moving_average(dN_d[:,i,j],n=nnn)
  rhod_avg[:,i,j] = moving_average(rho_d[:,i,j],n=nnn)
  ud_avg[:,i,j]   = moving_average(u_d[:,i,j],n=nnn)
  rd_avg[:,i,j]   = moving_average(rdust_d[:,i,j],n=nnn)
  lwhr_avg[:,i,j] = moving_average(lw_htrt_d[:,i,j],n=nnn)
  swhr_avg[:,i,j] = moving_average(sw_htrt_d[:,i,j],n=nnn)

## Plot settings (MUST CHANGE FROM MONTH TO MONTH)
######################################################################################
# Which Ls do you want to focus on?
Ls_ee= 153.9
Ls_e = 155.
l_1 = np.where(Ls - Ls_ee > 0.001)[0][0]
l_2 = np.where(Ls - Ls_e > 0.001)[0][0]
Ls = Ls[l_1:l_2]
n = l_2 - l_1

c = np.matrix('154. 0')#('0.5 45; 103.1 45; 242 2.5') # Dust storm mid-points [Ls Lat]

# Save destination

fpath = "/home/physastro/aes442/results/Dustruns/m%i/" % (amth)

#########################################################################################

######## TES dust files
# Zonal averaging
tau_d_z = d_d.sum(axis=2)/d_d.shape[2]

# Time averaging
nnnn=2
tau_d_avg=np.zeros((tau_d_z.shape[0]-(nnnn-1),tau_d_z.shape[1]))
for i in xrange(0,d_lat_s.shape[0]):
 tau_d_avg[:,i] = moving_average(tau_d_z[:,i],nnnn)

# first and last sols 
sol_a = int(np.round(669*(Ls_ee/360.)))
sol_s = int(np.round(669*(Ls_e/360.)))

tau_d_avg = tau_d_avg[sol_a:sol_s,:]
d_Ls_avg = np.linspace(Ls_ee,Ls_e,tau_d_avg.shape[0])
#########

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

## TEMP/WIND/TOPG map

# DATA
day = 1
hr = 96 # this is actually the tstep (t=96 is storm start)
lvl = 0

# variable[day][hour, elevation, lat, lon]
ut = ua[day][hr,lvl,:,:] - ub[day][hr,lvl,:,:]
vt = va[day][hr,lvl,:,:] - vb[day][hr,lvl,:,:]

#data = tempa[day][hr,lvl,:,:] - tempb[day][hr,lvl,:,:]
data = tsurfa[day][hr,:,:] - tsurfb[day][hr,:,:]
data2= presa[day][hr,:,:] - presb[day][hr,:,:]

# Longitude
major_ticksx = np.arange(np.floor(lon_t[0]), np.ceil(lon_t[-1]), 30)
minor_ticksx = np.arange(np.floor(lon_t[0]), np.ceil(lon_t[-1]), 10)

# Latitude
major_ticksy = np.arange(np.floor(lat_t[-1]), np.ceil(lat_t[0]), 30)
minor_ticksy = np.arange(np.floor(lat_t[-1]), np.ceil(lat_t[0]), 10)

## PLOT temperature/winds/topography
f, axarr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,10), dpi=100)
x = lon_t
y = lat_t

xlabel = 'Longitude / degrees'
ylabel = 'Latitude / degrees'
cblabel= 'Temperature difference / K'
plt.xlabel(xlabel, fontsize=14, labelpad=10)
plt.ylabel(ylabel, fontsize=14, labelpad=10)

# Main plot
ax = axarr.pcolormesh(x, y, data, cmap='RdBu_r', norm=MidPointNorm(midpoint=0.))

# Secondary plot
ax2 = axarr.quiver(x, y, ut, vt, scale=2**2, units='y', width=0.1)
aq = axarr.quiverkey(ax2, 0.815, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')

# Topography
lvls = [-5,0,5,10,15]
ax3 = axarr.contour(mola[1], mola[0], mola[2], lvls, colors='k')
 
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
plt.savefig("%stemp_uvwind_mola_latvslon.png" % (fpath), bbox_inches='tight')
plt.close('all')

## Temperature PLOT
temp_t = tsurf_avg_t.T
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
cb = f.colorbar(ax1, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
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
 
ax1 = axarr[0,0].pcolormesh(x, y, fslw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
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

ax3 = axarr[1,0].pcolormesh(x, y, fssw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
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
s_l = [-2.05, -6.12, 242.7]   # landing site marking on plot (actually for 244.7, Ls is messed up)
ticky_latlon = [60,10,30,10]  # tick settings [xmajor,xminor,ymajor,yminor] ticks
ticky_latalt = [60,10,20,10]
int_Ls = int(np.ceil(Ls.shape[0]/(12*Months)))

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
plt_timeseries(dustq_diff[l_1:,:,:], lon_t, lat_t, Ls_m[1][l_1:], 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg / kg', 3, '%sDustqdiff_latlon_tseries_ds1.png' % (fpath), mola)

alt_t = alt # Height of 20.9km
latll = 26
dustq_diff_altlon = dustqa[1][l_1:,:,latll,:] - dustqb[1][l_1:,:,latll,:]
temp_diff_altlon = tempa[1][l_1:,:,latll,:] - tempb[1][l_1:,:,latll,:]

plt_timeseries(temp_diff_altlon, lon_t, alt_t, Ls, 4,4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Temperature difference / K', int_Ls, '%stemp_altlon_tseries_ds1.png' % (fpath))
a
plt_timeseries(dustq_diff_altlon, lon_t, alt_t, Ls_m[1][l_1:], 4, 4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Dust MMR difference / kg / kg', 3, '%sdustq_altlon_tseries_ds1.png' % (fpath))

plt.close('all')

## IMPACT CALCULATIONS

# Target area 
llat1, llat2  = -22.5, 22.5
llon1, llon2  = -20., 20.
lalt1, lalt2  = 0., 8.

# Target time window
ts1, ts2      = 96, 108

lat_1, lat_2 = np.where(lat - llat2 >= 0.001)[0][-1]+1, np.where(lat - llat1 >= 0.001)[0][-1]+1
lon_1, lon_2 = np.where(lon - llon1 >= 0.001)[0][0]-1, np.where(lon - llon2 >= 0.001)[0][0]-1
alt_1, alt_2 = np.where(alt - lalt1 >= 0.001)[0][0], np.where(alt - lalt2 >= 0.001)[0][0]

alt_1 = 0

# Loop to compute impact
re_err = {}
re = {}
day = 1 

var_da = [tempa[1], presa[1], ua[1], va[1], rhoa[1], dustqa[1], dustNa[1], fluxsurflwa[1], fluxsurfswa[1], fluxtoplwa[1], fluxtopswa[1]]
var_db = [tempb[1], presb[1], ub[1], vb[1], rhob[1], dustqb[1], dustNb[1], fluxsurflwb[1], fluxsurfswb[1], fluxtoplwb[1], fluxtopswb[1]]

re[day] = np.zeros([len(var_da), (ts2-ts1)+1])
re_err[day] = np.zeros(len(var_da)) 
 
for n in xrange(0, len(var_da)):
 data_a = var_da[n]
 data_b = var_db[n]

 if len(data_a.shape)==4:

  m=0
  for j in xrange(ts1, ts2+1):
   aa = data_a[j,alt_1:alt_2,lat_1:lat_2,lon_1:lon_2].flatten() - data_b[j,alt_1:alt_2,lat_1:lat_2,lon_1:lon_2].flatten()
   a_ref = data_b[j,alt_1:alt_2,lat_1:lat_2,lon_1:lon_2].flatten()
   re[day][n,m] = np.linalg.norm(aa) / np.linalg.norm(a_ref)
   m=m+1
   
 else:
 
  m=0
  for j in xrange(ts1, ts2+1):
   aa = data_a[j,lat_1:lat_2,lon_1:lon_2].flatten() - data_b[j,lat_1:lat_2,lon_1:lon_2].flatten()
   a_ref = data_b[j,lat_1:lat_2,lon_1:lon_2].flatten()
   re[day][n,m] = np.linalg.norm(aa) / np.linalg.norm(a_ref)
   m=m+1
   
 re[day][(np.isnan(re[day])==True)] = 0.
 re_err[day][n] = sum(re[day][n,:]) / re[day][n,:].shape[0]

