# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors
import math

from mars_time import MarsTime
from scipy.io import *
from matplotlib import cm,ticker
from plt_timeseries import plt_timeseries, MidpointNormalize

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
ml = netcdf.netcdf_file('/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/mgcm-v6/datafile/surface.nc','r')

mola = {}

mola[0] = ml.variables['latitude'][:]
mola[1] = ml.variables['longitude'][:]
mola[2] = ml.variables['zMOL'][:]

# Number of months in comparison (always add 1 because of Python indexing)
Months = 13

# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,Months,1):
 mgcm = "MGCM_v5-1"
 rundira = "a_ds4"
 rundirb = "a_ref2"
 month = "m%i" % (i)
 filename = "diagfi.nc"
 
 a = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/%s/runs/%s/%s/%s" % (mgcm,rundira,month,filename),'r')
 b = netcdf.netcdf_file("/home/physastro/aes442/padata/alpha/users/aes442/mars_gcms/%s/runs/%s/%s/%s" % (mgcm,rundirb,month,filename),'r')

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
 taua[i] = a.variables['taudustvis'][:]
 rdusta[i] = a.variables['reffdust'][:]
 lw_htrta[i] = a.variables['lw_htrt'][:]
 sw_htrta[i] = a.variables['sw_htrt'][:]
 
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

#wd = 12 # writediagfi interval in runscript
#t = np.arange(1,669,float(1)/wd) # Create new time dimension variable 

## Ls vector
Ls_s = (Months-1)*30 # Number of solar longitudes for time vector for comparison
Ls = np.zeros((n))

# Method 1 creates a new Ls vector (it's neat, but doesn't line up with model times exactly)
#Ls = np.linspace(0,Ls_s,n)

# Method 2 grabs Ls's from model directly (has bugs, but can be ironed out)
p=0
for i in xrange(1,len(Ls_m)+1,1):
 gg = Ls_m[i]
 for j in xrange(gg.shape[0]):
  Ls[p] = gg[j]
  p = p + 1
Ls = np.roll(Ls,4)
Ls[-1] = Ls_s
Ls[:5] = [0,0.02,0.04,0.06,0.08]
print Ls.shape, Ls[:15], Ls[-15:]

## Create all other variables, with altitude dimension removed
ps_a, ps_b                 = np.zeros((n,lat.shape[0],lon.shape[0])), np.zeros((n,lat.shape[0],lon.shape[0]))
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
 
 # Dummy variables - reduce 4D data to 3D data by selecting surface only (time,lat,lon) 
 dmm2, dmm2b = dm2[:,1,:,:], dm2b[:,1,:,:]
 dmm3, dmm3b = dm3[:,1,:,:], dm3b[:,1,:,:]
 dmm5, dmm5b = dm5[:,1,:,:], dm5b[:,1,:,:]
 dmm6, dmm6b = dm6[:,1,:,:], dm6b[:,1,:,:]
 dmm7, dmm7b = dm7[:,1,:,:], dm7b[:,1,:,:] 
 dmm13, dmm13b = dm13[:,1,:,:], dm13b[:,1,:,:]
 dmm14, dmm14b = dm14[:,1,:,:], dm14b[:,1,:,:]
 dmm15, dmm15b = dm15[:,1,:,:], dm15b[:,1,:,:]

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
  
  m = m+1

tempaa, tempbb = {}, {}
uaa, ubb = {}, {}
dustqaa, dustqbb = {}, {}
dustNaa, dustNbb = {}, {}
rhoaa, rhobb = {}, {}
rdustaa, rdustbb = {}, {}
lw_htrtaa, lw_htrtbb = {}, {}
sw_htrtaa, sw_htrtbb = {}, {}

for k in xrange(1,len(psa)+1,1):
 dm2, dm2b = tempa[k], tempb[k]
 dm3, dm3b = ua[k], ub[k]
 dm5, dm5b = dustqa[k], dustqb[k]
 dm6, dm6b = dustNa[k], dustNb[k]
 dm7, dm7b = rhoa[k], rhob[k]
 dm13, dm13b = rdusta[k], rdustb[k]
 dm14, dm14b = lw_htrta[k], lw_htrtb[k]
 dm15, dm15b = sw_htrta[k], sw_htrtb[k]

 d2, d2b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d3, d3b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d5, d5b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d6, d6b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d7, d7b   = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d13, d13b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d14, d14b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))
 d15, d15b = np.zeros((len(tempa[k]),sigma.shape[0],lat.shape[0])), np.zeros((len(tempb[k]),sigma.shape[0],lat.shape[0]))

 for j in xrange(0,lon.shape[0],1):
  d2, d2b = d2 + dm2[:,:,:,j], d2b + dm2b[:,:,:,j]
  d3, d3b = d3 + dm3[:,:,:,j], d3b + dm3b[:,:,:,j]
  d5, d5b = d5 + dm5[:,:,:,j], d5b + dm5b[:,:,:,j]
  d6, d6b = d6 + dm6[:,:,:,j], d6b + dm6b[:,:,:,j]
  d7, d7b = d7 + dm7[:,:,:,j], d7b + dm7b[:,:,:,j]
  d13, d13b = d13 + dm13[:,:,:,j], d13b + dm13b[:,:,:,j]
  d14, d14b = d14 + dm14[:,:,:,j], d14b + dm14b[:,:,:,j]
  d15, d15b = d15 + dm15[:,:,:,j], d15b + dm15b[:,:,:,j]
  
 tempaa[k], tempbb[k]       = d2, d2b
 uaa[k], ubb[k]             = d3, d3b
 dustqaa[k], dustqbb[k]     = d5, d5b
 dustNaa[k], dustNbb[k]     = d6, d6b
 rhoaa[k], rhobb[k]         = d7, d7b
 rdustaa[k], rdustbb[k]     = d13, d13b
 lw_htrtaa[k], lw_htrtbb[k] = d14, d14b
 sw_htrtaa[k], sw_htrtbb[k] = d15, d15b
 
# Create variables now without longitude
temp_aa, temp_bb       = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
u_aa, u_bb             = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustq_aa, dustq_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
dustN_aa, dustN_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rho_aa, rho_bb         = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
rdust_aa, rdust_bb     = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
lw_htrt_aa, lw_htrt_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))
sw_htrt_aa, sw_htrt_bb = np.zeros((n,sigma.shape[0],lat.shape[0])), np.zeros((n,sigma.shape[0],lat.shape[0]))

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
 
 for i in xrange(dd2.shape[0]):
  temp_aa[m,:,:], temp_bb[m,:,:]   = dd2[i,:,:], dd2b[i,:,:]
  u_aa[m,:,:], u_bb[m,:,:]         = dd3[i,:,:], dd3b[i,:,:]
  dustq_aa[m,:,:], dustq_bb[m,:,:] = dd5[i,:,:], dd5b[i,:,:]
  dustN_aa[m,:,:], dustN_bb[m,:,:] = dd6[i,:,:], dd6b[i,:,:]
  rho_aa[m,:,:], rho_bb[m,:,:]     = dd7[i,:,:], dd7b[i,:,:]
  rdust_aa[m,:,:], rdust_bb[m,:,:] = dd13[i,:,:], dd13b[i,:,:]
  lw_htrt_aa[m,:,:], lw_htrt_bb[m,:,:] = dd14[i,:,:], dd14b[i,:,:]
  sw_htrt_aa[m,:,:], sw_htrt_bb[m,:,:] = dd15[i,:,:], dd15b[i,:,:]
  
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

# Initialise variables to calculate differences between reference run and dust storm run
dustq_diff   = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
dustN_diff   = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
temp_diff    = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
ps_diff      = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
rho_diff     = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
u_diff       = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
rdust_diff   = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
lw_htrt_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
sw_htrt_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
 
fluxsurf_lw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxsurf_sw_diff = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_lw_diff  = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))
fluxtop_sw_diff  = np.zeros((t.shape[0],lat.shape[0],lon.shape[0]))

dq_d      = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
dN_d      = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
t_d       = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
rho_d     = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
u_d       = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
rdust_d   = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
lw_htrt_d = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))
sw_htrt_d = np.zeros((t.shape[0],sigma.shape[0],lat.shape[0]))

# Calculate differences
for i in xrange(t.shape[0]):
 dustq_diff[i,:,:] = dustq_a[i,:,:] - dustq_b[i,:,:]
 dustN_diff[i,:,:] = dustN_a[i,:,:] - dustN_b[i,:,:]
 temp_diff[i,:,:] = temp_a[i,:,:] - temp_b[i,:,:]
 ps_diff[i,:,:] = ps_a[i,:,:] - ps_b[i,:,:]
 rho_diff[i,:,:] = rho_a[i,:,:] - rho_b[i,:,:]
 u_diff[i,:,:] = u_a[i,:,:] - u_b[i,:,:]
 rdust_diff[i,:,:] = rdust_a[i,:,:] - rdust_b[i,:,:]
 lw_htrt_diff[i,:,:] = lw_htrt_a[i,:,:] - lw_htrt_b[i,:,:]
 sw_htrt_diff[i,:,:] = sw_htrt_a[i,:,:] - sw_htrt_b[i,:,:]
  
 fluxsurf_lw_diff[i,:,:] = fluxsurf_lwa[i,:,:] - fluxsurf_lwb[i,:,:]
 fluxsurf_sw_diff[i,:,:] = fluxsurf_swa[i,:,:] - fluxsurf_swb[i,:,:]
 fluxtop_lw_diff[i,:,:] = fluxtop_lwa[i,:,:] - fluxtop_lwb[i,:,:]
 fluxtop_sw_diff[i,:,:] = fluxtop_swa[i,:,:] - fluxtop_swb[i,:,:]
 
 t_d[i,:,:] = temp_aa[i,:,:] - temp_bb[i,:,:]
 dq_d[i,:,:] = dustq_aa[i,:,:] - dustq_bb[i,:,:]
 dN_d[i,:,:] = dustN_aa[i,:,:] - dustN_bb[i,:,:]
 rho_d[i,:,:] = rho_aa[i,:,:] - rho_bb[i,:,:]
 u_d[i,:,:] = u_aa[i,:,:] - u_bb[i,:,:]
 rdust_d[i,:,:] = rdust_aa[i,:,:] - rdust_bb[i,:,:]
 lw_htrt_d[i,:,:] = lw_htrt_aa[i,:,:] - lw_htrt_bb[i,:,:]
 sw_htrt_d[i,:,:] = sw_htrt_aa[i,:,:] - sw_htrt_bb[i,:,:]

# Zonal averaging (time,lat)
temp_avg = np.zeros((t.shape[0],lat.shape[0]))
ps_avg = np.zeros((t.shape[0],lat.shape[0]))
u_avg = np.zeros((t.shape[0],lat.shape[0]))
rho_avg = np.zeros((t.shape[0],lat.shape[0]))
fssw_avg = np.zeros((t.shape[0],lat.shape[0]))
fslw_avg = np.zeros((t.shape[0],lat.shape[0]))
ftsw_avg = np.zeros((t.shape[0],lat.shape[0]))
ftlw_avg = np.zeros((t.shape[0],lat.shape[0]))
rdust_avg = np.zeros((t.shape[0],lat.shape[0]))
lw_htrt_avg = np.zeros((t.shape[0],lat.shape[0]))
sw_htrt_avg = np.zeros((t.shape[0],lat.shape[0]))
tau_a_avg, tau_b_avg = np.zeros((n,lat.shape[0])), np.zeros((n,lat.shape[0]))

for i in xrange(0,lon.shape[0],1):
 temp_avg = temp_avg + temp_diff[:,:,i]
 ps_avg = ps_avg + ps_diff[:,:,i]
 u_avg = u_avg + u_diff[:,:,i]
 rho_avg = rho_avg + rho_diff[:,:,i]
 fssw_avg = fssw_avg + fluxsurf_sw_diff[:,:,i]
 fslw_avg = fslw_avg + fluxsurf_lw_diff[:,:,i]
 ftsw_avg = ftsw_avg + fluxtop_sw_diff[:,:,i]
 ftlw_avg = ftlw_avg + fluxtop_lw_diff[:,:,i]
 tau_a_avg = tau_a_avg + tau_a[:,:,i]
 tau_b_avg = tau_b_avg + tau_b[:,:,i]
 rdust_avg = rdust_avg + rdust_diff[:,:,i]
 lw_htrt_avg = lw_htrt_avg + lw_htrt_diff[:,:,i]
 sw_htrt_avg = sw_htrt_avg + sw_htrt_diff[:,:,i]
 
temp_avg    = temp_avg/len(lon)
ps_avg      = ps_avg/len(lon)
u_avg       = u_avg/len(lon)
rho_avg     = rho_avg/len(lon)
fssw_avg    = fssw_avg/len(lon)
fslw_avg    = fslw_avg/len(lon)
ftsw_avg    = ftsw_avg/len(lon)
ftlw_avg    = ftlw_avg/len(lon)
tau_a_avg   = tau_a_avg/len(lon)
tau_b_avg   = tau_b_avg/len(lon)
rdust_avg   = rdust_avg/len(lon)
lw_htrt_avg = lw_htrt_avg/len(lon)
sw_htrt_avg = sw_htrt_avg/len(lon)

# Time moving-point average of Zonal average
nn=Ls.shape[0]/(Months-1) # Number of points to average over (30Ls, is about 1 month)
t_avg = Ls[:-(nn-1)]

temp_avg_t    = np.zeros((t_avg.shape[0],lat.shape[0]))
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

for i in xrange(0,lat.shape[0],1): 
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

# Calculate approximate HEIGHT from sigma (km)
alt = np.zeros((sigma.shape[0]))
for i in xrange(len(sigma)):
 alt[i] = -10.8*math.log(sigma[i])

c = np.matrix('1 45; 103 45; 243 2.5') # Dust storm mid-points [Ls Lat]

## PLOTS

# Common settings (ticks)
t_t = np.linspace(0,360,t_avg.shape[0])
t_tau = np.linspace(0,360,n)
lat_t = np.linspace(90,-90,lat.shape[0])
lon_t = np.linspace(-180,180,lon.shape[0])

# Solar longitude
major_ticksx = np.arange(0, Ls_s+1, 40)                                              
minor_ticksx = np.arange(0, Ls_s+1, 10)
# Latitude
major_ticksy = np.arange(-90, 91, 30)                                              
minor_ticksy = np.arange(-90, 91, 10)

## tau_ref and tau_ds PLOT
tau_ds = np.matrix.transpose(tau_a_avg)
tau_ref = np.matrix.transpose(tau_b_avg)

f, axarr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t_tau
y = lat_t
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'

cb_label = 'Dust optical depth / SI'
cb_label2 = 'Dust optical depth / SI'
 
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')

ax1 = axarr[0].pcolormesh(x, y, tau_ds, cmap='PuRd')
axarr[0].axis('tight') 
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 f.gca().add_artist(cmarker)
 axarr[0].plot(c[i,0],c[i,1],'o',color='y')
axarr[0].set_xticks(major_ticksx)
axarr[0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0].set_yticks(major_ticksy)                                                       
axarr[0].set_yticks(minor_ticksy, minor=True)
axarr[0].set_title('(a) Dust storm run', fontsize=18)
axarr[0].tick_params(axis='both', labelsize=14)

ax2 = axarr[1].pcolormesh(x, y, tau_ref, cmap='PuRd')
axarr[1].set_title('(b) Reference run', fontsize=18)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.54, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax1, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=14)                    # colorbar label

cbar_ax2 = f.add_axes([0.85, 0.1, 0.04, 0.36])  # [h_placement, v_placement, h_size, v_size]
cb2 = f.colorbar(ax2, cax=cbar_ax2, format='%.1f', extend='both') # double-edged colorbar
cb2.set_label('%s' % (cb_label2), fontsize=14)                    # colorbar label

plt.savefig("/home/physastro/aes442/CDOD_latvsLs_dsrunvsrefrun.png") 

## Temperature PLOT
temp_t = np.matrix.transpose(temp_avg_t)

fig = plt.figure(figsize=(10,10), dpi=100)
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(t_t,lat_t,temp_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.xlabel('Solar longitude / degrees',fontsize=18)
plt.ylabel('Latitude / degrees',fontsize=18)

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
ax.tick_params(axis='both', labelsize=14)

# Colour bar
cb = plt.colorbar(format='%.0f')
cb.set_label('Temperature difference / K')
tick_locator = ticker.MaxNLocator(nbins=16)
cb.locator = tick_locator

plt.axis('tight')
plt.savefig("/home/physastro/aes442/SurfTempDiff_LatvsTime_FY_uavg_tavg.png")

## Surface pressure and Atmospheric density at surface PLOT
ps_t = np.matrix.transpose(ps_avg_t)
rho_t = np.matrix.transpose(rho_avg_t)

f, axarr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
x = t_t
y = lat_t
xlabel = 'Solar longitude / degrees'
ylabel = 'Latitude / degrees'

cb_label = 'Pressure difference / Pa'
cb_label2 = 'Density difference / kg / m^-3'
 
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

# Zonal wind PLOT
u_t = np.matrix.transpose(u_avg_t)

fig = plt.figure( figsize=(10,10), dpi=100)
ax = fig.add_subplot(1,1,1)
plt.pcolormesh(t_t,lat_t,u_t,norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r')
plt.axis('tight')
plt.xlabel('Solar longitude / degrees',fontsize=16)
plt.ylabel('Latitude / degrees',fontsize=16)

for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 ax.plot(c[i,0],c[i,1],'o',color='y') 

ax.set_xticks(major_ticksx)
ax.set_xticks(minor_ticksx, minor=True)
ax.set_yticks(major_ticksy)
ax.set_yticks(minor_ticksy, minor=True)
 
cb = plt.colorbar(format='%.1f')
cb.set_label('Zonal wind velocity difference / ms^-1', rotation=270, labelpad=20)
tick_locator = ticker.MaxNLocator(nbins=7)
cb.locator = tick_locator
cb.update_ticks()

plt.savefig("/home/physastro/aes442/ZonalWindDiff_LatvsTime_FY_uavg_tavg.png")

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
cb_label = 'Surface radiative flux difference / W m^-2'     
# Common axis labels
f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')
 
ax1 = axarr[0,0].pcolormesh(x, y, fslw_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r', vmin=np.min(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)), vmax=np.max(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)))
axarr[0,0].axis('tight')
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 axarr[0,0].plot(c[i,0],c[i,1],'o',color='y')
axarr[0,0].set_xticks(major_ticksx)
axarr[0,0].set_xticks(minor_ticksx, minor=True)                                           
axarr[0,0].set_yticks(major_ticksy)                                                       
axarr[0,0].set_yticks(minor_ticksy, minor=True)
axarr[0,0].set_title('Incoming LW (IR) radiation (a)', fontsize=10)
axarr[0,0].tick_params(axis='both', labelsize=14)

ax2 = axarr[0,1].pcolormesh(x, y, ftlw_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r', vmin=np.min(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)), vmax=np.max(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)))
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 axarr[0,1].plot(c[i,0],c[i,1],'o',color='y')
axarr[0,1].set_title('Outgoing LW (IR) radiation (b)', fontsize=10)
axarr[0,1].tick_params(axis='both', labelsize=14)

ax3 = axarr[1,0].pcolormesh(x, y, fssw_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r', vmin=np.min(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)), vmax=np.max(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)))
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 axarr[1,0].plot(c[i,0],c[i,1],'o',color='y')
axarr[1,0].set_title('Incoming SW (VL) radiation (c)', fontsize=10)
axarr[1,0].tick_params(axis='both', labelsize=14)

ax4 = axarr[1,1].pcolormesh(x, y, ftsw_t, norm=MidpointNormalize(midpoint=0.), cmap='RdBu_r', vmin=np.min(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)), vmax=np.max(np.concatenate((fslw_t, ftlw_t, fssw_t, ftsw_t),axis=0)))
for i in xrange(len(c)):
 cmarker = plt.Circle((c[i,0],c[i,1],lat),0.2, color='black', fill=False)
 fig.gca().add_artist(cmarker)
 axarr[1,1].plot(c[i,0],c[i,1],'o',color='y')
axarr[1,1].set_title('Outgoing SW (VL) radiation (d)', fontsize=10)
axarr[1,1].tick_params(axis='both', labelsize=14)

# Colorbar creation and placement
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85, 0.1, 0.05, 0.8])  # [h_placement, v_placement, h_size, v_size]
cb = f.colorbar(ax3, cax=cbar_ax, format='%.1f', extend='both') # double-edged colorbar
cb.set_label('%s' % (cb_label), fontsize=14)                    # colorbar label

plt.savefig("/home/physastro/aes442/fluxes_latvsLs_zonavg_tavg.png")

### Time series plots
# settings
s_l = [-2.05, -6.12, 256.5]                      # landing site marking on plot (actually for 244.7, Ls is messed up)
ticky_latlon = [7,7,10]                          # tick settings [nbins,major_ticks,minor_ticks] (for lat/lon plots)
ticky_latalt = [7,7,10]                          # tick settings [nbins,major_ticks,minor_ticks] (for lat/alt plots)
int_Ls = int(math.ceil(Ls.shape[0]/(30*Months))) # interval of 1 Ls
ds2_Ls = (Ls.shape[0]/(30*(Months-1)))*(120.5)   # Start of DS2
ds3_Ls = (Ls.shape[0]/(30*(Months-1)))*(256)     # Start of DS3

# Dust particle size contours
rd_ds1, rd_ds2, rd_ds3 = {}, {}, {}
rd_ds1[0], rd_ds2[0], rd_ds3[0] = alt, alt, alt
rd_ds1[1], rd_ds2[1], rd_ds3[1] = lat_t, lat_t, lat_t
rd_ds1[2], rd_ds2[2], rd_ds3[2] = rd_avg, rd_avg[ds2_Ls:,:,:], rd_avg[ds3_Ls:,:,:]

# DS2 dust mmr average difference contours
dqd_ds = {}
dqd_ds[0] = alt
dqd_ds[1] = lat_t
dqd_ds[2] = dqd_avg[ds2_Ls:,:,:]

## Dust storm 1 Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff, lon_t, lat_t, Ls, 5,4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg/kg', int_Ls, 'Dustq_zontimeavg_tseries_ds1.png', mola)

#plt_timeseries(u_a, lon_t, lat_t, Ls, 5, 4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Zonal wind velocity / m/s', 36, 'uds_zontimeavg_tseries_ds1.png', mola)

#plt_timeseries(u_b, lon_t, lat_t, Ls, 5, 4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Zonal wind velocity / m/s', 36, 'uref_zontimeavg_tseries_ds1.png', mola)

## Dust storm 2 Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff[ds2_Ls:,:,:], lon_t, lat_t, Ls[ds2_Ls:], 5, 4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg/kg', int_Ls, 'Dustq_zontimeavg_tseries_ds2.png', mola)

## Dust storm 3 - Exo Mars  (-2.05145, -6.12423), Ls: 244.7, ToD: 13:15.
## Time series dustq (mmr) (time, lat, lon)
plt_timeseries(dustq_diff[ds3_Ls:,:,:], lon_t, lat_t, Ls[ds3_Ls:], 5, 4, ticky_latlon, 'Longitude / degrees', 'Latitude / degrees', 'Ls: ', 'Dust MMR difference / kg/kg', int_Ls, 'Dustq_zontimeavg_tseries_exomars.png', mola, s_l)

## DS1 temp timeseries (time, alt, lat)
#plt_timeseries(td_avg, lat_t, alt, t_av, 3, 3, ticky_latalt, 'Latitude / degrees', 'Altitude above surface / km', 'Ls: ', 'Temperature difference / K', int_Ls, 'tempdiff_zontimeavg_tseries_ds1.png',rd_ds1)
#plt_timeseries(dqd_avg, lat_t, alt, t_av, 5, 4, ticky_latalt, 'Latitude / degrees', 'Altitude above surface / km', 'Ls: ', 'Dust MMR difference / kg/kg', int_Ls, 'dustdiff_zontimeavg_tseries_ds1.png')

## DS2 temp timeseries (time, alt, lat)
#plt_timeseries(td_avg[ds2_Ls:,:,:], lat_t, alt, t_av[ds2_Ls:], 3, 3, ticky_latalt, 'Latitude / degrees', 'Altitude above surface / km', 'Ls: ', 'Temperature difference / K', int_Ls, 'tempdiff_zontimeavg_tseries_ds2.png', rd_ds2)
plt_timeseries(dqd_avg[ds2_Ls:,:,:], lat_t, alt, t_av[ds2_Ls:], 5, 4, ticky_latalt, 'Latitude / degrees', 'Altitude above surface / km', 'Ls: ', 'Dust MMR difference / kg/kg', int_Ls, 'dustdiff_zontimeavg_tseries_ds2.png', rd_ds2)

## DS3 temp timeseries (time, alt, lat)
#plt_timeseries(td_avg[ds3_Ls:,:,:], lat_t, alt, t_av[ds3_Ls:], 3, 3, ticky_latalt, 'Latitude / degrees', 'Altitude above surface / km', 'Ls: ', 'Temperature difference / K', int_Ls, 'tempdiff_zontimeavg_tseries_exomars.png', rd_ds3)
#plt_timeseries(dqd_avg[ds3_Ls:,:,:], lat_t, alt, t_av[ds3_Ls:], 5, 4, ticky_latalt, 'Latitude / degrees', 'Altitude above surface / km', 'Ls: ', 'Dust MMR difference / kg/kg', int_Ls, 'dustdiff_zontimeavg_tseries_exomars.png')

## DS 3 daily time series temp/heatrate (time, alt)
#plt_timeseries(lwhr_avg[ds2_Ls:,:,:], lat_t, alt, t_av[ds2_Ls:], 3, 3, ticky_latalt, 'Latitude / degrees', 'Altitude above surface / km', 'Ls: ', 'Long wave heating rate difference / K / sol^-1', int_Ls, 'lwhr_dqd_zontimeavg_tseries_ds2.png',dqd_ds)

#fix trigger to change from contour to pcolormesh
#fix trigger to label plots properly
#fix trigger to change the magnitude of units in colorbar
