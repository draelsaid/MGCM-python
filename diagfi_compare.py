# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016

import matplotlib as mpl
mpl.use('Agg') # removes need for X-Server (sshing graphics in linux). For qsub only.

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

from mars_time import MarsTime
from scipy.io import *
from matplotlib import cm,ticker
#from plt_timeseries import *
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Use tex font
#rc('text',usetex=True)
# Change all fonts to 'Computer Modern'
#rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})

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
ua, ub = {}, {}
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

# Grab topography from surface.nc or mola32.nc file
ml = netcdf.netcdf_file('/padata/alpha/users/aes442/mars_gcms/mgcm-v6/datafile/surface.nc','r')

mola = {}

mola[0] = ml.variables['latitude'][:]
mola[1] = ml.variables['longitude'][:]
mola[2] = ml.variables['zMOL'][:]

# Import data from Luca's TES dust files for comparison
a = netcdf.netcdf_file('/padata/alpha/users/aes442/mars_gcms/mgcm-v6/datafile/dust_MY28.nc','r')

d_lat_s = a.variables['latitude'][:]
d_lon_s = a.variables['longitude'][:]
d_t   = a.variables['Time'][:]
d_d   = a.variables['dustop'][:]

d_lat = np.linspace(-90,90,d_lat_s.shape[0])
d_lon = np.linspace(-180,180,d_lon_s.shape[0])

# Number of months in comparison (always add 1 because of Python indexing)
Months = 13

# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,Months):
 mgcm = "MGCM_v5-1"
 rundira = "a_ds7"
 rundirb = "a_ref4"
 month = "m%i" % (i)
 filename = "diagfi.nc"
 
 a = netcdf.netcdf_file("/padata/alpha/users/aes442/mars_gcms/%s/runs/%s/%s/%s" % (mgcm,rundira,month,filename),'r')
 b = netcdf.netcdf_file("/padata/alpha/users/aes442/mars_gcms/%s/runs/%s/%s/%s" % (mgcm,rundirb,month,filename),'r')

 lat = a.variables['lat'][:]
 lon = a.variables['lon'][:]
 sigma = a.variables['sigma'][:]
 t_m = a.variables['time'][:]
 Ls_m[i] = a.variables['Ls'][:]
 
 psa[i] = a.variables['ps'][:]
 presa[i]= a.variables['pressure'][:]
 tempa[i] = a.variables['temp'][:]
 ua[i] = a.variables['u'][:]
 dustqa[i] = a.variables['dustq'][:]
 dustNa[i] = a.variables['dustN'][:]
 rhoa[i] = a.variables['rho'][:]
 fluxsurflwa[i] = a.variables['fluxsurf_lw'][:]
 fluxsurfswa[i] = a.variables['fluxsurf_sw'][:]
 fluxtoplwa[i] = a.variables['fluxtop_lw'][:]
 fluxtopswa[i] = a.variables['fluxtop_sw'][:]
 taua[i] = a.variables['taudustvis'][:]
 rdusta[i] = a.variables['reffdust'][:]
 lw_htrta[i] = a.variables['lw_htrt'][:]
 sw_htrta[i] = a.variables['sw_htrt'][:]
 
 psb[i] = b.variables['ps'][:]
 presb[i]= b.variables['pressure'][:]
 tempb[i] = b.variables['temp'][:]
 ub[i] = b.variables['u'][:]
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
for i in xrange(1,len(psa)+1,1): # len(psa) gives the number of months
 n = n + len(dustqa[i])          # len(dustqa[i]) gives the number of time steps in each month. Different variable used as a cross-check of dimension consistency.
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

## Create all other variables, with altitude dimension removed
ps_a, ps_b                 = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
pres_a, pres_b             = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
temp_a, temp_b             = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
u_a, u_b                   = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustq_a, dustq_b           = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
dustN_a, dustN_b           = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
rho_a, rho_b               = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_lwa, fluxsurf_lwb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxsurf_swa, fluxsurf_swb = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_lwa, fluxtop_lwb   = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
fluxtop_swa, fluxtop_swb   = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
tau_a, tau_b               = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
rdust_a, rdust_b           = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
lw_htrt_a, lw_htrt_b       = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
sw_htrt_a, sw_htrt_b       = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))

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
 dm12, dm12b = taua[k], taub[k]
 dm13, dm13b = rdusta[k], rdustb[k]
 dm14, dm14b = lw_htrta[k], lw_htrtb[k]
 dm15, dm15b = sw_htrta[k], sw_htrtb[k]
 dm16, dm16b = presa[k], presb[k]
 
 # Dummy variables - reduce 4D data to 3D data by selecting surface only (time,lat,lon) 
 dmm2, dmm2b = dm2[:,1,:,:], dm2b[:,1,:,:]
 dmm3, dmm3b = dm3[:,1,:,:], dm3b[:,1,:,:]
 dmm5, dmm5b = dm5[:,1,:,:], dm5b[:,1,:,:]
 dmm6, dmm6b = dm6[:,1,:,:], dm6b[:,1,:,:]
 dmm7, dmm7b = dm7[:,1,:,:], dm7b[:,1,:,:] 
 dmm13, dmm13b = dm13[:,1,:,:], dm13b[:,1,:,:]
 dmm14, dmm14b = dm14[:,1,:,:], dm14b[:,1,:,:]
 dmm15, dmm15b = dm15[:,1,:,:], dm15b[:,1,:,:]
 dmm16, dmm16b = dm16[:,1,:,:], dm16b[:,1,:,:]

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
  tau_a[m,:,:], tau_b[m,:,:]              = dm12[i,:,:], dm12b[i,:,:]
  rdust_a[m,:,:], rdust_b[m,:,:]          = dmm13[i,:,:], dmm13b[i,:,:]
  lw_htrt_a[m,:,:], lw_htrt_b[m,:,:]      = dmm14[i,:,:], dmm14b[i,:,:]
  sw_htrt_a[m,:,:], sw_htrt_b[m,:,:]      = dmm15[i,:,:], dmm15b[i,:,:]
  pres_a[m,:,:], pres_b[m,:,:]            = dmm16[i,:,:], dmm16b[i,:,:]
  
  m = m+1

tempaa, tempbb = {}, {}
uaa, ubb = {}, {}
dustqaa, dustqbb = {}, {}
dustNaa, dustNbb = {}, {}
rhoaa, rhobb = {}, {}
rdustaa, rdustbb = {}, {}
lw_htrtaa, lw_htrtbb = {}, {}
sw_htrtaa, sw_htrtbb = {}, {}
presaa, presbb = {}, {}

for k in xrange(1,len(psa)+1,1):
 dm2, dm2b = tempa[k], tempb[k]
 dm3, dm3b = ua[k], ub[k]
 dm5, dm5b = dustqa[k], dustqb[k]
 dm6, dm6b = dustNa[k], dustNb[k]
 dm7, dm7b = rhoa[k], rhob[k]
 dm13, dm13b = rdusta[k], rdustb[k]
 dm14, dm14b = lw_htrta[k], lw_htrtb[k]
 dm15, dm15b = sw_htrta[k], sw_htrtb[k]
 dm16, dm16b = presa[k], presb[k]

 d2, d2b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d3, d3b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d5, d5b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d6, d6b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d7, d7b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d13, d13b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d14, d14b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d15, d15b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d16, d16b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))

 for j in xrange(0,lon.shape[0],1):
  d2, d2b = d2 + dm2[:,:,:,j], d2b + dm2b[:,:,:,j]
  d3, d3b = d3 + dm3[:,:,:,j], d3b + dm3b[:,:,:,j]
  d5, d5b = d5 + dm5[:,:,:,j], d5b + dm5b[:,:,:,j]
  d6, d6b = d6 + dm6[:,:,:,j], d6b + dm6b[:,:,:,j]
  d7, d7b = d7 + dm7[:,:,:,j], d7b + dm7b[:,:,:,j]
  d13, d13b = d13 + dm13[:,:,:,j], d13b + dm13b[:,:,:,j]
  d14, d14b = d14 + dm14[:,:,:,j], d14b + dm14b[:,:,:,j]
  d15, d15b = d15 + dm15[:,:,:,j], d15b + dm15b[:,:,:,j]
  d16, d16b = d16 + dm16[:,:,:,j], d16b + dm16b[:,:,:,j]
  
 tempaa[k], tempbb[k]       = d2, d2b
 uaa[k], ubb[k]             = d3, d3b
 dustqaa[k], dustqbb[k]     = d5, d5b
 dustNaa[k], dustNbb[k]     = d6, d6b
 rhoaa[k], rhobb[k]         = d7, d7b
 rdustaa[k], rdustbb[k]     = d13, d13b
 lw_htrtaa[k], lw_htrtbb[k] = d14, d14b
 sw_htrtaa[k], sw_htrtbb[k] = d15, d15b
 presaa[k], presbb[k] = d15, d15b
 
# Create variables now without longitude
temp_aa, temp_bb       = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
u_aa, u_bb             = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustq_aa, dustq_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustN_aa, dustN_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rho_aa, rho_bb         = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rdust_aa, rdust_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
lw_htrt_aa, lw_htrt_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
sw_htrt_aa, sw_htrt_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
pres_aa, pres_bb       = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))

m=0
for k in xrange(1,len(psa)+1,1):
 dd2, dd2b = tempaa[k], tempbb[k]
 dd3, dd3b = uaa[k], ubb[k]
 dd5, dd5b = dustqaa[k], dustqbb[k]
 dd6, dd6b = dustNaa[k], dustNbb[k]
 dd7, dd7b = rhoaa[k], rhobb[k]
 dd13, dd13b = rdustaa[k], rdustbb[k]
 dd14, dd14b = lw_htrtaa[k], lw_htrtbb[k]
 dd15, dd15b = sw_htrtaa[k], sw_htrtbb[k]
 dd16, dd16b = presaa[k], presbb[k]
 
 for i in xrange(dd2.shape[0]):
  temp_aa[m,:,:], temp_bb[m,:,:]   = dd2[i,:,:], dd2b[i,:,:]
  u_aa[m,:,:], u_bb[m,:,:]         = dd3[i,:,:], dd3b[i,:,:]
  dustq_aa[m,:,:], dustq_bb[m,:,:] = dd5[i,:,:], dd5b[i,:,:]
  dustN_aa[m,:,:], dustN_bb[m,:,:] = dd6[i,:,:], dd6b[i,:,:]
  rho_aa[m,:,:], rho_bb[m,:,:]     = dd7[i,:,:], dd7b[i,:,:]
  rdust_aa[m,:,:], rdust_bb[m,:,:] = dd13[i,:,:], dd13b[i,:,:]
  lw_htrt_aa[m,:,:], lw_htrt_bb[m,:,:] = dd14[i,:,:], dd14b[i,:,:]
  sw_htrt_aa[m,:,:], sw_htrt_bb[m,:,:] = dd15[i,:,:], dd15b[i,:,:]
  pres_aa[m,:,:], pres_bb[m,:,:] = dd16[i,:,:], dd16b[i,:,:]
  
  m = m+1

# Zonal averaging (time,alt,lat)
temp_aa, temp_bb   = temp_aa/len(lon), temp_bb/len(lon)
u_aa, u_bb         = u_aa/len(lon), u_bb/len(lon) 
dustq_aa, dustq_bb = dustq_aa/len(lon), dustq_bb/len(lon)
dustN_aa, dustN_bb = dustN_aa/len(lon), dustN_bb/len(lon) 
rho_aa, rho_bb     = rho_aa/len(lon), rho_bb/len(lon)
rdust_aa, rdust_bb = rdust_aa/len(lon), rdust_bb/len(lon)
lw_htrt_aa, lw_htrt_bb = lw_htrt_aa/len(lon), lw_htrt_bb/len(lon)
sw_htrt_aa, sw_htrt_bb = sw_htrt_aa/len(lon), sw_htrt_bb/len(lon)
pres_aa, pres_bb = pres_aa/len(lon), pres_bb/len(lon)

# Calculate differences
dustq_diff = dustq_a - dustq_b
dustN_diff = dustN_a - dustN_b
temp_diff = temp_a - temp_b
ps_diff = ps_a - ps_b
rho_diff = rho_a - rho_b
u_diff = u_a - u_b
rdust_diff = rdust_a - rdust_b
lw_htrt_diff = lw_htrt_a - lw_htrt_b
sw_htrt_diff = sw_htrt_a - sw_htrt_b
pres_diff = pres_a - pres_b
  
fluxsurf_lw_diff = fluxsurf_lwa - fluxsurf_lwb
fluxsurf_sw_diff = fluxsurf_swa - fluxsurf_swb
fluxtop_lw_diff = fluxtop_lwa - fluxtop_lwb
fluxtop_sw_diff = fluxtop_swa - fluxtop_swb

t_d = temp_aa - temp_bb
dq_d = dustq_aa - dustq_bb
dN_d = dustN_aa - dustN_bb
rho_d = rho_aa - rho_bb
u_d = u_aa - u_bb
rdust_d = rdust_aa - rdust_bb
lw_htrt_d = lw_htrt_aa - lw_htrt_bb
sw_htrt_d = sw_htrt_aa - sw_htrt_bb
pres_d = pres_aa - pres_bb

# Zonal averaging (time,lat)

temp_avg = np.sum(temp_a,axis=2)/temp_a.shape[2] - np.sum(temp_b,axis=2)/temp_b.shape[2]
pres_avg = np.sum(pres_a,axis=2)/pres_a.shape[2] - np.sum(pres_b,axis=2)/pres_b.shape[2]
ps_avg = np.sum(ps_a,axis=2)/ps_a.shape[2] - np.sum(ps_b,axis=2)/ps_b.shape[2]
u_avg = np.sum(u_a,axis=2)/u_a.shape[2] - np.sum(u_b,axis=2)/u_b.shape[2]
rho_avg = np.sum(rho_a,axis=2)/rho_a.shape[2] - np.sum(rho_b,axis=2)/rho_b.shape[2]
fssw_avg = np.sum(fluxsurf_swa,axis=2)/fluxsurf_swa.shape[2] - np.sum(fluxsurf_swb,axis=2)/fluxsurf_swb.shape[2]
fslw_avg = np.sum(fluxsurf_lwa,axis=2)/fluxsurf_lwa.shape[2] - np.sum(fluxsurf_lwb,axis=2)/fluxsurf_lwb.shape[2]
ftsw_avg = np.sum(fluxtop_swa,axis=2)/fluxtop_swa.shape[2] - np.sum(fluxtop_swb,axis=2)/fluxtop_swb.shape[2]
ftlw_avg = np.sum(fluxtop_lwa,axis=2)/fluxtop_lwa.shape[2] - np.sum(fluxtop_lwb,axis=2)/fluxtop_lwb.shape[2]
tau_a_avg = np.sum(tau_a,axis=2)/tau_a.shape[2]
tau_b_avg = np.sum(tau_b,axis=2)/tau_b.shape[2]
rdust_avg = np.sum(rdust_a,axis=2)/rdust_a.shape[2] - np.sum(rdust_b,axis=2)/rdust_b.shape[2]
lw_htrt_avg = np.sum(lw_htrt_a,axis=2)/lw_htrt_a.shape[2] - np.sum(lw_htrt_b,axis=2)/lw_htrt_b.shape[2]
sw_htrt_avg = np.sum(sw_htrt_a,axis=2)/sw_htrt_a.shape[2] - np.sum(sw_htrt_b,axis=2)/sw_htrt_b.shape[2]

# Time moving-point average of Zonal average
nn=Ls.shape[0]/((Months-1)/3) # Number of points to average over
t_avg = Ls[:-(nn-1)]

temp_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
ps_avg_t      = np.zeros((t_avg.shape[0],lat.shape[0]))
pres_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
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
 pres_avg_t[:,i]    = moving_average(pres_avg[:,i],n=nn)

############ NOW NEED TIME AVERAGE FOR NEW DATA ###################
nnn=nn
t_av = Ls[:-(nnn-1)]

td_avg    = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
dqd_avg   = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
dNd_avg   = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
rhod_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
ud_avg    = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
rd_avg    = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
lwhr_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
swhr_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))
presd_avg  = np.zeros((t_av.shape[0],sigma.shape[0],lat.shape[0]))

for j in xrange(0,lat.shape[0],1):
 for i in xrange(0,sigma.shape[0],1):
  td_avg[:,i,j]   = moving_average(t_d[:,i,j],n=nnn)
  dqd_avg[:,i,j]  = moving_average(dq_d[:,i,j],n=nnn)
  dNd_avg[:,i,j]  = moving_average(dN_d[:,i,j],n=nnn)
  rhod_avg[:,i,j] = moving_average(rho_d[:,i,j],n=nnn)
  ud_avg[:,i,j]   = moving_average(u_d[:,i,j],n=nnn)
  rd_avg[:,i,j]   = moving_average(rdust_d[:,i,j],n=nnn)
  lwhr_avg[:,i,j] = moving_average(lw_htrt_d[:,i,j],n=nnn)
  swhr_avg[:,i,j] = moving_average(sw_htrt_d[:,i,j],n=nnn)
  presd_avg[:,i,j] = moving_average(pres_d[:,i,j],n=nnn)

# Calculate approximate HEIGHT from sigma (km)
alt = np.zeros((sigma.shape[0]))
for i in xrange(len(sigma)):
 alt[i] = -10.8*np.log(sigma[i])

c = np.matrix('4.01 45; 103.1 45; 242 2.5') # Dust storm mid-points [Ls Lat]

fpath = "/padata/alpha/users/aes442/results/Dustruns/" % (amth)

print "PLOTTING....."
## TES dust files
# Zonal averaging
nnnn=20
tau_d_z = d_d.sum(axis=2)/d_d.shape[2] 

tau_d_avg=np.zeros((tau_d_z.shape[0]-(nnnn-1),tau_d_z.shape[1]))
for i in xrange(0,d_lat_s.shape[0]):
 tau_d_avg[:,i] = moving_average(tau_d_z[:,i],nnnn)
d_Ls_avg = np.linspace(0,Ls_s,d_t.shape[0]-(nnnn-1))

## PLOTS
# Common settings (ticks)
t_t = np.linspace(0,Ls_s,t_avg.shape[0])
t_tau = np.linspace(0,Ls_s,n)
lat_t = np.linspace(90,-90,lat.shape[0])
lon_t = np.linspace(-180,180,lon.shape[0])

# Solar longitude
major_ticksx = np.arange(0, Ls_s+1, 40)                                              
minor_ticksx = np.arange(0, Ls_s+1, 10)
# Latitude
major_ticksy = np.arange(-90, 91, 30)                                              
minor_ticksy = np.arange(-90, 91, 10)

## tau_ref, tau_ds, tau_tes PLOT
tau_ds = np.matrix.transpose(tau_a_avg)
tau_ref = np.matrix.transpose(tau_b_avg)
tau_TES = np.matrix.transpose(tau_d_avg)

f, axarr = plt.subplots(3, 1, figsize=(12,12), dpi=100)
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

ax1 = axarr[0].pcolormesh(x, y, tau_ds, cmap='gist_rainbow_r',vmax=2.3,vmin=0)
axarr[0].axis('tight') 
for i in xrange(len(c)):
 axarr[0].plot(c[i,0],c[i,1],'o',color='y',markersize=10)
axarr[0].set_xticks(major_ticksx)
axarr[0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0].set_yticks(major_ticksy)                                                       
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('(a) Dust storm run', fontsize=14)
axarr[0].tick_params(axis='both', labelsize=14)

ax2 = axarr[1].pcolormesh(x, y, tau_ref, cmap='gist_rainbow_r',vmax=3,vmin=0)
axarr[1].axis('tight')
axarr[1].set_title('(b) Reference run', fontsize=14)

ax3 = axarr[2].pcolormesh(xx,yy,tau_TES, cmap='gist_rainbow_r',vmax=3,vmin=0)
axarr[2].axis('tight')
axarr[2].set_title('(c) Mars year 28', fontsize=14)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8])    # [h_placement, v_placement, h_size, v_size]
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

plt.savefig("%sCDOD_latvsLs_dsrunvsrefrun.png" % (fpath)) 

## Temperature PLOT
temp_t = np.matrix.transpose(temp_avg_t)

fig = plt.figure(figsize=(10,10), dpi=100)
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(t_t,lat_t,temp_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
plt.xlabel('Solar longitude / degrees',fontsize=18)
plt.ylabel('Latitude / degrees',fontsize=18)

# Extra Markers
for i in xrange(len(c)):
 ax.plot(c[i,0],c[i,1],'o',color='y',markersize=10)
 
# Ticks
ax.set_xticks(major_ticksx)                                                       
ax.set_xticks(minor_ticksx, minor=True)                                       
ax.set_yticks(major_ticksy)                                                       
ax.set_yticks(minor_ticksy, minor=True)
ax.tick_params(axis='both', labelsize=14)

# Colour bar
cb = plt.colorbar(format='%.1f')
cb.set_label('Temperature difference / K')
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator

plt.axis('tight')
plt.savefig("%sSurfTempDiff_LatvsTime_FY_uavg_tavg.png" % (fpath))

## Atmospheric pressure and density at surface PLOT
ps_t = np.matrix.transpose(pres_avg_t)
rho_t = np.matrix.transpose(rho_avg_t)

f, axarr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t_t
y = lat_t
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'

cb_label = 'Atmospheric pressure difference / Pa'
cb_label2 = 'Atmospheric density difference / kg / m^3'
 
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')

ax1 = axarr[0].pcolormesh(x, y, ps_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
axarr[0].axis('tight') 
for i in xrange(len(c)):
 axarr[0].plot(c[i,0],c[i,1],'o',color='y',markersize=10)
axarr[0].set_xticks(major_ticksx)
axarr[0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0].set_yticks(major_ticksy)                                                       
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('(a)', fontsize=18)
axarr[0].tick_params(axis='both', labelsize=14)

ax2 = axarr[1].pcolormesh(x, y, rho_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
for i in xrange(len(c)):
 axarr[1].plot(c[i,0],c[i,1],'o',color='y',markersize=10)
axarr[1].set_title('(b)', fontsize=18)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.54, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax1, cax=cbar_ax, format='%.0f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=14)                    # colorbar label

cbar_ax2 = f.add_axes([0.85, 0.1, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb2 = f.colorbar(ax2, cax=cbar_ax2, format='%.1e', extend='both') # double-edged colorbar
cb2.set_label('%s' % (cb_label2), fontsize=14)                    # colorbar label

plt.savefig("%sPresDensDiff_LatvsLs_zonavg_tavg.png" % (fpath)) 

# Zonal wind PLOT
u_t = np.matrix.transpose(u_avg_t)

fig = plt.figure( figsize=(10,10), dpi=100)
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(t_t,lat_t,u_t,norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Solar longitude / degrees',fontsize=16)
plt.ylabel('Latitude / degrees',fontsize=16)

for i in xrange(len(c)):
 ax.plot(c[i,0],c[i,1],'o',color='y',markersize=10) 
ax.set_xticks(major_ticksx)
ax.set_xticks(minor_ticksx, minor=True)
ax.set_yticks(major_ticksy)
ax.set_yticks(minor_ticksy, minor=True)
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Zonal wind velocity difference / m / s', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=7)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("%sZonalWindDiff_LatvsTime_FY_uavg_tavg.png" % (fpath))

# ALL FLUXES on one plot
fslw_t = np.matrix.transpose(fslw_avg_t) # Incoming long wave (IR) radiation
ftlw_t = np.matrix.transpose(ftlw_avg_t) # Outgoing long wave (IR) radiation
fssw_t = np.matrix.transpose(fssw_avg_t) # Incoming short wave (VL) radiation
ftsw_t = np.matrix.transpose(ftsw_avg_t) # Outgoing short wave (VL) radiation
    
f, axarr = plt.subplots(2, 2, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t_t
y = lat_t
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'
cb_label = 'Radiative flux difference / W / m^2' 
    
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')
 
ax1 = axarr[0,0].pcolormesh(x, y, fslw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')# vmin=np.min((np.min(fslw_t),np.min(ftlw_t),np.min(fssw_t),np.min(ftsw_t))), vmax=np.max((np.max(fslw_t),np.max(ftlw_t),np.max(fssw_t),np.max(ftsw_t))))
axarr[0,0].axis('tight')
for i in xrange(len(c)):
 axarr[0,0].plot(c[i,0],c[i,1],'o',color='y',markersize=10)
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
for i in xrange(len(c)):
 axarr[0,1].plot(c[i,0],c[i,1],'o',color='y',markersize=10)
axarr[0,1].set_title('Outgoing flux at top (LW) (b)', fontsize=10)
axarr[0,1].tick_params(axis='both', labelsize=10)
dv2 = make_axes_locatable(axarr[0,1])
cax2 = dv2.append_axes("right",size="5%",pad=0.05)
cb2 = f.colorbar(ax2,cax=cax2, format='%.1f', extend='both')
cb2.set_label('%s' % (cb_label), fontsize=10)

ax3 = axarr[1,0].pcolormesh(x, y, fssw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
for i in xrange(len(c)):
 axarr[1,0].plot(c[i,0],c[i,1],'o',color='y',markersize=10)
axarr[1,0].set_title('Incident flux at surface (SW) (c)', fontsize=10)
axarr[1,0].tick_params(axis='both', labelsize=10)
dv3 = make_axes_locatable(axarr[1,0])
cax3 = dv3.append_axes("right",size="5%",pad=0.05)
cb3 = f.colorbar(ax3,cax=cax3, format='%.1f', extend='both')
cb3.set_label('%s' % (cb_label), fontsize=10)

ax4 = axarr[1,1].pcolormesh(x, y, ftsw_t, norm=MidPointNorm(midpoint=0.), cmap='RdBu_r')
for i in xrange(len(c)):
 axarr[1,1].plot(c[i,0],c[i,1],'o',color='y',markersize=10)
axarr[1,1].set_title('Outgoing flux at top (SW) (d)', fontsize=10)
axarr[1,1].tick_params(axis='both', labelsize=10)
dv4 = make_axes_locatable(axarr[1,1])
cax4 = dv4.append_axes("right",size="5%",pad=0.05)
cb4 = f.colorbar(ax4,cax=cax4, format='%.1f', extend='both')
cb4.set_label('%s' % (cb_label), fontsize=10)   

plt.savefig("%sfluxes_latvsLs_zonavg_tavg.png" % (fpath))

### Time series plots
# settings
s_l = [-2.05, -6.12, 256.5]                      # landing site marking on plot (actually for 244.7, Ls is messed up)
ticky_latlon = [7,7,10]                          # tick settings [nbins,major_ticks,minor_ticks] (for lat/lon plots)
ticky_latalt = [7,7,10]                          # tick settings [nbins,major_ticks,minor_ticks] (for lat/alt plots)
int_Ls = int(np.ceil(Ls.shape[0]/(30*(Months-1)))) # Ls interval
int_Ls2 = int(np.ceil(Ls.shape[0]/(30*(Months-1)))/5) # Ls interval2
ds1_Ls = np.where(np.absolute(Ls-4)==np.min(np.absolute(Ls-4)))[0][0]
ds2_Ls = np.where(np.absolute(Ls-102)==np.min(np.absolute(Ls-102)))[0][0]   # Start of DS2
ds3_Ls = np.where(np.absolute(Ls-241)==np.min(np.absolute(Ls-241)))[0][0]   # Start of DS3
ds2_Lsa = np.where(np.absolute(t_av-102)==np.min(np.absolute(t_av-102)))[0][0] # Start of DS2 (t_av vect)
ds3_Lsa = np.where(np.absolute(t_av-242)==np.min(np.absolute(t_av-242)))[0][0]     # Start of DS3 (t_av vect)

# DS2 dust mmr average difference contours
dqd_ds = {}
dqd_ds[0] = alt
dqd_ds[1] = lat_t
dqd_ds[2] = dqd_avg[ds2_Lsa:,:,:]

exit()

## Dust storm 1 Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff[ds1_Ls:,:,:], lon_t, lat_t, Ls[ds1_Ls:], 4,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sDustq_zontimeavg_tseries_ds1.png' % (fpath), mola)

## Dust storm 2 Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff[ds2_Ls:,:,:], lon_t, lat_t, Ls[ds2_Ls:], 4, 4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sDustq_zontimeavg_tseries_ds2.png' % (fpath), mola)

## Dust storm 3 - Exo Mars  (-2.05145, -6.12423), Ls: 244.7, ToD: 13:15.
## Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff[ds3_Ls:,:,:], lon_t, lat_t, Ls[ds3_Ls:], 4, 4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sDustq_zontimeavg_tseries_ds3.png' % (fpath), mola)

alt_t = np.linspace(0,21,alt.shape[0]) # Height of 20.9km (approx)
dustq_diff_altlon = dustqa[1][ds1_Ls:,:,8,:] - dustqb[1][ds1_Ls:,:,8,:]
temp_diff_altlon = tempa[1][ds1_Ls:,:,8,:] - tempb[1][ds1_Ls:,:,8,:]

plt_timeseries(temp_diff_altlon, lon_t, alt_t, Ls[ds1_Ls:], 4,4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Temperature difference / K ', int_Ls, '%stemp_altlon_tseries_ds1.png' % (fpath))

plt_timeseries(dustq_diff_altlon, lon_t, alt_t, Ls[ds1_Ls:], 4,4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sdustq_altlon_tseries_ds1.png' % (fpath))

dustq_diff_altlon = dustqa[8][50:,:,17,:] - dustqb[8][50:,:,17,:]
temp_diff_altlon = tempa[8][50:,:,17,:] - tempb[8][50:,:,17,:]

plt_timeseries(temp_diff_altlon, lon_t, alt_t, Ls_m[9][27:], 4,4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Temperature difference / K ', int_Ls, '%stemp_altlon_tseries_ds3.png' % (fpath))

plt_timeseries(dustq_diff_altlon, lon_t, alt_t, Ls_m[9][27:], 4,4, ticky_latalt, 'Longitude / degrees', 'Altitude / km', 'Ls: ', 'Dust MMR difference / kg / kg', int_Ls, '%sdustq_altlon_tseries_ds3.png' % (fpath))

plt.close()
#fix trigger to change from contour to pcolormesh
#fix trigger to label plots properly
#fix trigger to change the magnitude of units in colorbar
