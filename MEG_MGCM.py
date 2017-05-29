# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016

import matplotlib as mpl
#mpl.use('Agg') # removes need for X-Server (sshing graphics in linux). For qsub only.

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as colors

from mars_time import MarsTime
from pylab import *
from scipy.io import *
from fmcd import call_mcd,julian
from matplotlib import pyplot

# Use tex font
#rc('text',usetex=True)
# Change all fonts to 'Computer Modern'
#rc('font',**{'family':'sans-serif','sans-serif':['Computer Modern']})

# Initialise dictionaries - due to data size
Ls_m = {}
psa, psb = {}, {}
presa, presb = {}, {} 
tempa, tempb = {}, {} 
vb, ub = {}, {}
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

# Number of months in comparison (always add 1 because of Python indexing)
Months = 13

# This loop assigns the data in both directories to variables here. This is done for each month. The result is a dictionary of dictionaries. One dictionary containing a dictionary for every month.
for i in xrange(1,Months):
 mgcm = "MGCM_v5-1"
 rundirb = "T31_9_Marsyears/MY31"
 month = "m%i" % (i)
 filename = "diagfi.nc"

 b = netcdf.netcdf_file("/padata/alpha/users/aes442/RUNS/R-%s/%s/%s/%s" % (mgcm,rundirb,month,filename),'r')

# b = netcdf.netcdf_file("/padata/mars/users/rmc429/testruns/v5dustlifting/1131_201604_v5-1hT31L25p24WDf/%s/%s" % (month, filename),'r')

 lat = b.variables['lat'][:]
 lon = b.variables['lon'][:]
 sigma = b.variables['sigma'][:]
 t_m = b.variables['time'][:]
# Ls_m[i] = b.variables['Ls'][:]

 psb[i] = b.variables['ps'][:]
# presb[i]= b.variables['pressure'][:]
 tempb[i] = b.variables['temp'][:]
 ub[i] = b.variables['u'][:]
 vb[i] = b.variables['v'][:]
 rhob[i] = b.variables['rho'][:]

# for a run without pressure
 presb[i] = tempb[i]*0.
 for k in xrange(sigma.shape[0]):
  presb[i][:,k,:,:] = psb[i][:,:,:]*sigma[k]
 
print ("Latitude: %i ||" % (lat.shape)), ("Longitude: %i ||" % (lon.shape)), ("Model levels: %i ||" % (sigma.shape))

# Get time dimension length
n = 0
for i in xrange(1,len(psb)+1,1): # len(psa) gives the number of months
 n = n + len(tempb[i])          # len(dustqa[i]) gives the number of time steps in each month. Different variable used as a cross-check of dimension consistency.
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

# Calculate approximate HEIGHT from sigma (km)
alt = np.zeros((sigma.shape[0]))
for i in xrange(len(sigma)):
 alt[i] = -10.8*np.log(sigma[i])

print "MCD PART................"

# Command abbreviations for sol and Ls conversions
ls_sol=MarsTime().ls_sol # Takes Ls and returns sol
sol_ls=MarsTime().sol_ls # Takes sol and returns Ls

## 1.3 Position (latitude, longitude) (deg) 
latt = 0.0
lont = 0.0

lat_m = np.arange(latt-5.,latt+5.+1.,1)
lon_m = np.arange(lont-5.,lont+5.+1.,1)

## 1.4 Dust and solar scenario
dust=np.array([1,7,8])

dust_sc = {}
for i in xrange(len(dust)):
 if (dust[i]==1):
  dust_sc[i] = "Climatology avg solar EUV"
 if (dust[i]==2): 
  dust_sc[i] = "Climatology min solar EUV"
 if (dust[i]==3):
  dust_sc[i] = "Climatology max solar EUV"
 if (dust[i]==4):
  dust_sc[i] = "Dust storm (tau=5) min EUV"
 if (dust[i]==5):
  dust_sc[i] = "Dust storm (tau=5) avg EUV"
 if (dust[i]==6):
  dust_sc[i] = "Dust storm (tau=5) max EUV"
 if (dust[i]==7):
  dust_sc[i] = "Dustier (warm) than clim. max EUV"
 if (dust[i]==8):
  dust_sc[i] = "Clearer (cold) than clim. min EUV"
 if (dust[i]==24) or (dust[i]==25) or (dust[i]==26) or (dust[i]==27) or (dust[i]==28) or (dust[i]==29) or (dust[i]==30) or (dust[i]==31):
  dust_sc[i] = ("MarsYear: %s" % (dust[i]))

dset='/padata/alpha/users/aes442/mcd/MCD_DATA/'
gwlength=16000
datekey=1
zkey=2
alt_m=concatenate((np.arange(0,15000,50),np.arange(15000,100100,100)),axis=0)
#alt_m = 1000.*np.array([  9.99499977e-01,   9.98098314e-01,   9.95562136e-01,
#         9.90976512e-01,   9.82737064e-01,   9.68096793e-01,
#         9.42590177e-01,   8.99628282e-01,   8.31171811e-01,
#         7.31004179e-01,   6.00909054e-01,   4.55260634e-01,
#         3.16687942e-01,   2.03942463e-01,   1.23567358e-01,
#         7.16656074e-02,   4.03222963e-02,   2.21809167e-02,
#         1.19544938e-02,   6.29041623e-03,   3.19916755e-03,
#         1.54171500e-03,   6.77837757e-04,   2.49743753e-04,
#         5.62741079e-05])

hrkey=1
varmodel=1 # Var model IS NOT REQUIRED for MEGT.
extvarkeys = np.ones(100)

# (1) Initialise Ls array
Ls_m = np.array([1.,360.])

# convert to sols, create vector, then convert each sol back to Ls again (for 1sol intervals, but in Ls)
sol_a = ls_sol(float(Ls_m[0]))
sol_b = ls_sol(float(Ls_m[-1]))

Ls_m = np.arange(sol_a,sol_b+1)
for i in xrange(len(Ls_m)):
 Ls_m[i] = sol_ls(Ls_m[i])
 
Ls_m[0] = 1.
Ls_m[-1]=360.

# (2) Initialise local time of day array
tod_c = np.array([0,24])
tb=int(tod_c[0])
te=int(tod_c[1])
xloct = np.arange(tb,te+1,3)

###################################################################
#                             MCD Call                            #
###################################################################

p_dim = len(Ls_m)*len(xloct)#*len(lat_m)*len(lon_m) # profile dimension

# Initialise arrays
mcd_temp    = np.zeros((p_dim,len(dust),len(alt_m)))
mcd_merwind = np.zeros((p_dim,len(dust),len(alt_m)))
mcd_zonwind = np.zeros((p_dim,len(dust),len(alt_m)))
mcd_dens    = np.zeros((p_dim,len(dust),len(alt_m)))
mcd_pres    = np.zeros((p_dim,len(dust),len(alt_m)))

print "Data being retrieved from the MCD is of size: ", mcd_temp.shape

# MCD loop
seedin=0.
for l in xrange(len(dust)):
 m=0
 for j in xrange(len(Ls_m)):
  print "Ls: ", Ls_m[j]
  for k in xrange(len(xloct)):
   for n in xrange(len(alt_m)):
     
    (pres, dens, temp, zonwind, merwind, \
    meanvar, extvar, seedout, ierr) \
    = \
    call_mcd(zkey,alt_m[n],lont,latt,hrkey, \
    datekey,Ls_m[j],xloct[k],dset,dust[l], \
    varmodel,seedin,gwlength,extvarkeys )
     
    mcd_merwind[m,l,n] = merwind
    mcd_zonwind[m,l,n] = zonwind
    mcd_dens[m,l,n] = dens
    mcd_temp[m,l,n] = temp
    mcd_pres[m,l,n] = pres
     
   m = m + 1      # Ls,xloct counter
   
#################################################################################################

mcd_temp[(mcd_temp == -999.) | (mcd_temp == 0.)] = np.NaN
mcd_dens[(mcd_dens == -999.) | (mcd_dens == 0.)] = np.NaN
mcd_zonwind[(mcd_zonwind == -999.) | (mcd_zonwind == 0.)] = np.NaN
mcd_temp[(mcd_merwind == -999.) | (mcd_merwind == 0.)] = np.NaN
mcd_pres[(mcd_pres == -999.) | (mcd_pres == 0.)] = np.NaN

alt_m = alt_m/1000.

print "MCD CALL OVER. PLOTTING....."

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

#data
lvl = 0
latt1 = 17
latt2 = 0
lonn = 36

y = alt
y2 = alt_m

# Common axis labels
cmap = mpl.cm.hsv

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,2):
 for k in xrange(50):
  if (k==1 and i==1):
   ax = axr.plot(ub[j][k,:,latt1,lonn], y, alpha=0.5, linewidth=1.5, color=cmap(0.5), label="MGCM - MY31")
  else:
   ax = axr.plot(ub[j][k,:,latt1,lonn], y, alpha=0.5, linewidth=1.5, color=cmap(0.5))
  
for k in xrange(len(dust)):
 for i in xrange(50):
  if (k==1 and i==1): 
   ax = axr.plot(mcd_zonwind[i,k,:], y2, alpha=0.3, linewidth=1.5, color=cmap(1.), label="MCD - Scenarios") 
  else:
   ax = axr.plot(mcd_zonwind[i,k,:], y2, alpha=0.3, linewidth=1.5, color=cmap(1.))
  
legend = plt.legend(loc='best', ncol=1, fontsize=9)
for l in legend.get_lines():
 l.set_alpha(1)
 
plt.axis([-400, 250, 0, 100])
axr.set_xlabel('Zonal wind velocity / m/s', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('u_profile.png')

print "done"
exit()

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(vb[j].shape[0]):
  ax = axr.plot(vb[j][k,:,latt1,lonn], y, alpha=0.15, linewidth=1.5, color=cmap(1))
  
plt.axis([-250, 250, 0, 100])
axr.set_xlabel('Meridional wind velocity / m/s', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('v_profile.png')

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(presb[j].shape[0]):
  ax = axr.plot(presb[j][k,:,latt1,lonn], y, alpha=0.15, linewidth=1.5, color=cmap(1))
  
plt.axis([0, 750, 0, 100])
axr.set_xlabel('Pressure / Pa', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('pressure_profile.png')

f,axr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)
for j in xrange(1,Months):
 for k in xrange(rhob[j].shape[0]):
  ax = axr.plot(rhob[j][k,:,latt1,lonn], y, alpha=0.15, linewidth=1.5, color=cmap(1))
  
plt.axis([0, 0.02, 0, 100])
axr.set_xlabel('Density / kg/$m^3$', fontsize=12)
axr.set_ylabel('Height above Mars areoid / km', fontsize=12)
plt.savefig('density_profile.png')

