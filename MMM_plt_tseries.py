import numpy as np
import matplotlib.pyplot as plt

from MidPointNorm import *

def MMM_plt_tseries(data, x, y, z, rows, cols, xlabel, ylabel, sp_titles, cb_label, zint, filename, d2_c=0):

    f, axarr = plt.subplots(rows, cols, sharex=True, sharey=True, figsize=(12,12), dpi=200)
    #plt.setp(axarr.flat, aspect=1.0, adjustable='box-forced')
    
# Assign labels
    xlabel = xlabel
    ylabel = ylabel
    cb_label = cb_label
    f.text(0.5, 0.04, '%s' % (xlabel), fontsize=18, ha='center')
    f.text(0.06, 0.5, '%s' % (ylabel), fontsize=18, va='center', rotation='vertical')
    
# Ticks
    major_ticksx = np.arange(np.floor(x[0]), np.ceil(x[-1]), 20)
    minor_ticksx = np.arange(np.floor(x[0]), np.ceil(x[-1]), 5)
    major_ticksy = np.arange(np.floor(y[0]), np.ceil(y[-1]), 20)
    minor_ticksy = np.arange(np.floor(y[0]), np.ceil(y[-1]), 5)

# Subplots
    m=0
    for i in xrange(rows):
     for j in xrange(cols):

      ax = axarr[i,j].pcolormesh(x, y, data[m*zint,:,:], cmap='RdBu_r', norm=MidPointNorm(midpoint=0.), vmin=np.min(data))
      axarr[i,j].axis('tight') 
      axarr[i,j].set_title('%s %0.1f' % (sp_titles, z[m*zint]), fontsize=10)
      axarr[i,j].set_xticks(major_ticksx)
      axarr[i,j].set_yticks(major_ticksy)   
      axarr[i,j].set_xticks(minor_ticksx, minor=True)                                      
      axarr[i,j].set_yticks(minor_ticksy, minor=True)
      axarr[i,j].tick_params(axis='both', labelsize=8, pad=10)
      
      if (d2_c!=0):
       lvls = [-5,0,5,10,15]
       ax2 = axarr[i,j].contour(d2_c[1], d2_c[0], d2_c[2], lvls, colors='k')
      
      m=m+1

# Colorbar creation and placement
    if (np.max(data)<=0.09):
     f.subplots_adjust(right=0.8)
     cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8]) # [h_plc, v_plc, h_size, v_size]
     cb = f.colorbar(ax, cax=cbar_ax, format='%.2e', extend='both')
     cb.set_label('%s' % (cb_label), fontsize=16) # colorbar label
    elif (np.max(data)<=100) & (np.min(data)>=0.1):
     f.subplots_adjust(right=0.8)
     cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8]) # [h_plc, v_plc, h_size, v_size]
     cb = f.colorbar(ax, cax=cbar_ax, format='%.1f', extend='both')
     cb.set_label('%s' % (cb_label), fontsize=16) # colorbar label
    else:
     f.subplots_adjust(right=0.8)
     cbar_ax = f.add_axes([0.85, 0.1, 0.04, 0.8]) # [h_plc, v_plc, h_size, v_size]
     cb = f.colorbar(ax, cax=cbar_ax, format='%.0f', extend='both')
     cb.set_label('%s' % (cb_label), fontsize=16) # colorbar label     
     
    plt.savefig("%s" % (filename), bbox_inches='tight', dpi=200)
