# this module is used for dealing with the time issue on Mars 
# and it is worth to expand this routine to cover other
# time-related stuffs. 
# sol_ls & ls_sol are based on the idl routines 
#
# This python module is written by LJ Ruan
# usage: 
# from mars_time import MarsTime
# sol_Ls=MarsTime().sol_ls

import math

class MarsTime:
    def __init__(self):
        self.year_day=668.6
	self.peri_day=485.0
	self.timeperi=1.905637
	self.e_elips=0.093358
	self.pi=3.1415927
	self.degrad=57.295779

    def sol_ls(self,sol):
	zz=(sol-self.peri_day)/self.year_day
	zanom=2.*self.pi*(zz-float(int(zz+0.5)))
	xref=abs(zanom)
	#The equation zx0 - e * sin (zx0) = xref, solved by Newton
	zx0=xref+self.e_elips*math.sin(xref)
	for iter in range(10):
  	    zdx=-(zx0-self.e_elips*math.sin(zx0)-xref)/(1.-self.e_elips*math.cos(zx0))
            if abs(zdx) <= 1.e-7:
	       break
  	    zx0=zx0+zdx
        else:
            print 'solution did not converge in 10 times'
        zx0=zx0+zdx
        if zanom < 0.:
           zx0=-zx0
	zteta=2.*math.atan(math.sqrt((1.+self.e_elips)/(1.-self.e_elips))*math.tan(zx0/2.))
	ls=zteta-self.timeperi
        if ls < 0.:
           ls=ls+2.*self.pi
	if ls >= 2.*self.pi:
	   ls=ls-2.*self.pi
	ls=self.degrad*ls

	return ls

    def ls_sol(self, ls):
	if (abs(ls) < 1.0e-5):
           if (ls >= 0.0):
              return 0.0
	   else:
	      return self.year_day

        zteta=ls/self.degrad+self.timeperi
	zx0=2.0*math.atan(math.tan(0.5*zteta)/math.sqrt((1.+self.e_elips)/(1.-self.e_elips)))
	xref=zx0-self.e_elips*math.sin(zx0)
	zz=xref/(2.*self.pi)
	sol=zz*self.year_day+self.peri_day
        if (sol < 0.0):
           sol=sol+self.year_day
        elif (sol >= self.year_day):
           sol=sol-self.year_day

 	return sol


