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

Months = 2   # No. of months
amth = 1     # Actual month 

for i in xrange(1,Months):
 mgcm = "MGCM_v5-1"
 rundira = "a_ds8"
 rundirb = "a_ref4"
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
n = len(dustqa[1])
print ("Total time steps: %i" % (n))
# these are the middle of the relevant gridboxes
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

f,axr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)

lvl = 20

temp_ds = tempa[1][:,:lvl,8:11,9]
temp_ref = tempb[1][:,:lvl,8:11,9]

y = alt[:lvl]

# ticks
tmajor_ticksx = np.arange(140, 240+1, 20)                                              
tminor_ticksx = np.arange(140, 240+1, 5)

tmajor_ticksy = np.arange(0, max(y)+10, 10)                                              
tminor_ticksy = np.arange(0, max(y)+10, 5)

# Common axis labels
ylabel = 'Height above Mars areoid / km'
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=14, va='center', rotation='vertical')

# Colors
cmap1 = mpl.cm.PuOr
cmap2 = mpl.cm.PuBu_r

t_days = 4
t_strt = 89

########################
for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(temp_ref.shape[2]):
  if (92<=k<=98) | (104<=k<=110) | (116<=k<=122) | (128<=k<=134)| (128<=k<=134) | (140<=k<=146) | (152<=k<=158) | (164<=k<=170) | (176<=k<=182) | (188<=k<=194):
   ax2 = axr[0].plot(temp_ref[k,:,j],y, color=cmap1(40), alpha=0.6, label='Day' if i==0 else "")
  else:
   ax2 = axr[0].plot(temp_ref[k,:,j],y,color=cmap2(40), alpha=0.45, label='Night' if i==0 else "")
axr[0].set_title('Reference run', fontsize=12)
########################
for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(temp_ds.shape[2]):
  if (92<=k<=98) | (104<=k<=110) | (116<=k<=122) | (128<=k<=134)| (128<=k<=134) | (140<=k<=146) | (152<=k<=158) | (164<=k<=170) | (176<=k<=182) | (188<=k<=194):
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

plt.axis([140, 240, 0, np.max(y)])
plt.suptitle('Sols 7-10. Lat: 40N-50N. Lon: 135W.')
plt.savefig('temp_profile.png')

#######################

f,axr = plt.subplots(2, 1, sharex=True, sharey=True, figsize=(10,12), dpi=200)

lvl = 20

u_ds = ua[1][:,:lvl,8:11,9]
u_ref = ub[1][:,:lvl,8:11,9]

y = alt[:lvl]

# ticks
umajor_ticksx = np.arange(-40, 140+1, 20)                                              
uminor_ticksx = np.arange(-40, 140+1, 5)

umajor_ticksy = np.arange(0, max(y)+10, 10)                                              
uminor_ticksy = np.arange(0, max(y)+10, 5)

# Common axis labels
ylabel = 'Height above Mars areoid / km'
f.text(0.06, 0.5, '%s' % (ylabel), fontsize=14, va='center', rotation='vertical')

for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(u_ref.shape[2]):
  if (92<=k<=98) | (104<=k<=110) | (116<=k<=122) | (128<=k<=134)| (128<=k<=134) | (140<=k<=146) | (152<=k<=158) | (164<=k<=170) | (176<=k<=182) | (188<=k<=194):
   ax4 = axr[0].plot(u_ref[k,:,j],y, color=cmap1(40), alpha=0.6, label='Day' if i==0 else "")
  else:
   ax4 = axr[0].plot(u_ref[k,:,j],y,color=cmap2(40), alpha=0.45, label='Night' if i==0 else "")
axr[0].set_title('Reference run', fontsize=12)
####################
for k in xrange(t_strt,t_strt+12*t_days+1):
 for j in xrange(u_ds.shape[2]):
  if (92<=k<=98) | (104<=k<=110) | (116<=k<=122) | (128<=k<=134) | (140<=k<=146) | (152<=k<=158) | (164<=k<=170) | (176<=k<=182) | (188<=k<=194):
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

plt.axis([-40, 140, 0, np.max(y)])
axr[1].set_title('Dust storm run', fontsize=12)
axr[1].set_xlabel('Zonal wind velocity / m/s', fontsize=12)
plt.suptitle('Sols 7-10. Lat: 40N-50N. Lon: 135W.')
plt.savefig('u_profile.png')



