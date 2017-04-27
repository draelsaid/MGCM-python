# Compares NetCDF data from the Mars GCM for Full Mars Year by combining monthly output of diagfi.nc files
# Adam El-Said 08/2016
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pylab as pltt
import matplotlib as mpl

from pylab import *
from scipy.io import *

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

Months = 12   # No. of months
amth = 30     # Actual month 

for i in xrange(1,Months):
 mgcm = "MGCM_v5-1"
 rundirb = "a_ref4"
 month = ("m%s" % (i)) # CHANGE
 filename = "diagfi.nc"
 
 a = netcdf.netcdf_file("/padata/alpha/users/aes442/RUNS/R-%s/%s/%s/%s" % (mgcm,rundira,month,filename),'r')
 
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
 rhoa[i]        = a.variables['rho'][:]
 fluxsurflwa[i] = a.variables['fluxsurf_lw'][:]
 fluxsurfswa[i] = a.variables['fluxsurf_sw'][:]
 fluxtoplwa[i]  = a.variables['fluxtop_lw'][:]
 fluxtopswa[i]  = a.variables['fluxtop_sw'][:]
 taua[i]        = a.variables['taudustvis'][:]
 rdusta[i]      = a.variables['reffdust'][:]

# Calculate approximate HEIGHT from sigma (km)
alt = np.zeros((sigma.shape[0]))
for i in xrange(len(sigma)):
 alt[i] = -10.8*np.log(sigma[i])

print "Latitude: %i || Longitude: %i || Model levels: %i => Alt Min:%.3f | Alt Max:%.3f | Alt half: %.3f " % (lat.shape[0],lon.shape[0],sigma.shape[0],alt[0],alt[-1],alt[18])
alt_half=18 # 47.8km

# Get time dimension length
n = len(dustqa[1])
print ("Total time steps: %i" % (n))

f,axr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)

lvl = 20
latt1 = 18
latt2 = 19
lonn = 36

temp_ds = tempa[1][:,:lvl,latt1:latt2,lonn]
temp_ref = tempb[1][:,:lvl,latt1:latt2,lonn]

y = alt[:lvl]

tmax = 270
tmin = 120

# ticks
tmajor_ticksx = np.arange(tmin, tmax+1, 20)                                              
tminor_ticksx = np.arange(tmin, tmax+1, 5)

tmajor_ticksy = np.arange(0, max(y)+10, 10)                                              
tminor_ticksy = np.arange(0, max(y)+10, 5)

# Common axis labels
ylabel = 'Height above Mars areoid / km'
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=14, va='center', rotation='vertical')

# Colors
cmap1 = mpl.cm.PuOr
cmap2 = mpl.cm.PuBu_r

t_days = 1
sol_start = 84        # sol start (from midnight)
t_strt = sol_start+3  # daytime starts

########################
for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(temp_ref.shape[2]):
  if (t_strt<=k<=t_strt+6) | (t_strt+12<=k<=t_strt+18) | (t_strt+24<=k<=t_strt+30) | (t_strt+36<=k<=t_strt+40)| (t_strt+48<=k<=t_strt+54) | (t_strt+60<=k<=t_strt+66) | (t_strt+72<=k<=t_strt+78) | (t_strt+84<=k<=t_strt+90) | (t_strt+96<=k<=t_strt+102) | (t_strt+108<=k<=t_strt+114):
   ax2 = axr[0].plot(temp_ref[k,:,j],y, color=cmap1(40), alpha=0.6, label='Day' if i==0 else "")
  else:
   ax2 = axr[0].plot(temp_ref[k,:,j],y,color=cmap2(40), alpha=0.45, label='Night' if i==0 else "")
axr[0].set_title('Reference run', fontsize=12)
########################
for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(temp_ds.shape[2]):
  if (t_strt<=k<=t_strt+6) | (t_strt+12<=k<=t_strt+18) | (t_strt+24<=k<=t_strt+30) | (t_strt+36<=k<=t_strt+40)| (t_strt+48<=k<=t_strt+54) | (t_strt+60<=k<=t_strt+66) | (t_strt+72<=k<=t_strt+78) | (t_strt+84<=k<=t_strt+90) | (t_strt+96<=k<=t_strt+102) | (t_strt+108<=k<=t_strt+114):
   ax1 = axr[1].plot(temp_ds[k,:,j],y, color=cmap1(40), alpha=0.6, label='Day' if i==0 else "")
  else:
   ax1 = axr[1].plot(temp_ds[k,:,j],y,color=cmap2(40), alpha=0.45, label='Night' if i==0 else "")
axr[1].set_title('Dust storm run', fontsize=12)
axr[1].set_xlabel('Temperature / K', fontsize=12)
######################
axr[1].set_xticks(tmajor_ticksx)
axr[1].set_xticks(tminor_ticksx, minor=True)
axr[1].set_yticks(tmajor_ticksy)
axr[1].set_yticks(tminor_ticksy, minor=True)

axr[0].set_xticks(tmajor_ticksx)
axr[0].set_xticks(tminor_ticksx, minor=True)
axr[0].set_yticks(tmajor_ticksy)
axr[0].set_yticks(tminor_ticksy, minor=True)

plt.axis([tmin, tmax, 0, np.max(y)])
plt.suptitle('Sols 9-13, month 6. Lat: 10S-10N. Lon: 0.')
plt.savefig('temp_profile.png')

#######################

f,axr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)

u_ds = ua[1][:,:lvl,latt1:latt2,lonn]
u_ref = ub[1][:,:lvl,latt1:latt2,lonn]

y = alt[:lvl]

umax = 60
umin = -90

# ticks
umajor_ticksx = np.arange(umin, umax+1, 20)                                              
uminor_ticksx = np.arange(umin, umax+1, 5)

umajor_ticksy = np.arange(0, max(y)+10, 10)                                              
uminor_ticksy = np.arange(0, max(y)+10, 5)

# Common axis labels
ylabel = 'Height above Mars areoid / km'
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=14, va='center', rotation='vertical')

for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(u_ref.shape[2]):
  if (t_strt<=k<=t_strt+6) | (t_strt+12<=k<=t_strt+18) | (t_strt+24<=k<=t_strt+30) | (t_strt+36<=k<=t_strt+40)| (t_strt+48<=k<=t_strt+54) | (t_strt+60<=k<=t_strt+66) | (t_strt+72<=k<=t_strt+78) | (t_strt+84<=k<=t_strt+90) | (t_strt+96<=k<=t_strt+102) | (t_strt+108<=k<=t_strt+114):
   ax4 = axr[0].plot(u_ref[k,:,j],y, color=cmap1(40), alpha=0.6, label='Day' if i==0 else "")
  else:
   ax4 = axr[0].plot(u_ref[k,:,j],y,color=cmap2(40), alpha=0.45, label='Night' if i==0 else "")
axr[0].set_title('Reference run', fontsize=12)
####################
for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(u_ds.shape[2]):
  if (t_strt<=k<=t_strt+6) | (t_strt+12<=k<=t_strt+18) | (t_strt+24<=k<=t_strt+30) | (t_strt+36<=k<=t_strt+40)| (t_strt+48<=k<=t_strt+54) | (t_strt+60<=k<=t_strt+66) | (t_strt+72<=k<=t_strt+78) | (t_strt+84<=k<=t_strt+90) | (t_strt+96<=k<=t_strt+102) | (t_strt+108<=k<=t_strt+114):
   ax3 = axr[1].plot(u_ds[k,:,j],y, color=cmap1(40), alpha=0.6, label='Day' if i==0 else "")
  else:
   ax3 = axr[1].plot(u_ds[k,:,j],y,color=cmap2(40), alpha=0.45, label='Night' if i==0 else "")

axr[1].set_xticks(umajor_ticksx)
axr[1].set_xticks(uminor_ticksx, minor=True)
axr[1].set_yticks(umajor_ticksy)
axr[1].set_yticks(uminor_ticksy, minor=True)

axr[0].set_xticks(umajor_ticksx)
axr[0].set_xticks(uminor_ticksx, minor=True)
axr[0].set_yticks(umajor_ticksy)
axr[0].set_yticks(uminor_ticksy, minor=True)

plt.axis([umin, umax, 0, np.max(y)])
axr[1].set_title('Dust storm run', fontsize=12)
axr[1].set_xlabel('Zonal wind velocity / m/s', fontsize=12)
plt.suptitle('Sols 9-13, month 6. Lat: 10S-10N. Lon: 0.')
plt.savefig('u_profile.png')



