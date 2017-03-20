'''
========================================================
Demonstration of advanced quiver and quiverkey functions
========================================================

Known problem: the plot autoscaling does not take into account
the arrows, so those on the boundaries are often out of the picture.
This is *not* an easy problem to solve in a perfectly general way.
The workaround is to manually expand the axes.
'''
import matplotlib.pyplot as plt
import numpy as np
from numpy import ma

def plt_quiver(X, Y, U, V, xlabel, ylabel, fname):

	f, axarr = plt.subplots(1, 1, sharex=True, sharey=True, figsize=(12,12), dpi=100)
# Assign labels
	xlabel = xlabel
	ylabel = ylabel
	f.text(0.5, 0.04, '%s' % (xlabel), fontsize=12, ha='center')
	f.text(0.06, 0.5, '%s' % (ylabel), fontsize=12, va='center', rotation='vertical')
	
	X, Y = np.meshgrid(X, Y)
#	U = np.cos(X)
#	V = np.sin(Y)

	axarr.set_title('Arrows scale with plot width, not view', fontsize=12)
	Q = axarr.quiver(X, Y, U, V, units='width', linewidths=1.)
	q = axarr.quiverkey(Q, 0.9, 0.9, 1, r'$2 \frac{m}{s}$', labelpos='E', coordinates='figure')
        axarr.axis('tight') 
        
	plt.savefig("%s" % (fname), bbox_inches='tight')
	
#plt.figure()
#plt.title("pivot='mid'; every third arrow; units='inches'")
#Q = plt.quiver(X[::3, ::3], Y[::3, ::3], U[::3, ::3], V[::3, ::3],
#               pivot='mid', units='inches')
#qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E', coordinates='figure')
#plt.scatter(X[::3, ::3], Y[::3, ::3], color='r', s=5)

#plt.figure()
#plt.title("pivot='tip'; scales with x view")
#M = np.hypot(U, V)
#Q = plt.quiver(X, Y, U, V, M, units='x', pivot='tip', width=0.022, scale=1 / 0.15)
#qk = plt.quiverkey(Q, 0.9, 0.9, 1, r'$1 \frac{m}{s}$', labelpos='E',
#                   coordinates='figure')
#plt.scatter(X, Y, color='k', s=5)

	
