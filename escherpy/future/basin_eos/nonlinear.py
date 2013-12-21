#!/usr/bin/env python

'''
Non-linear fitting with Levenberg-Marquardt
to a Birch--Murnaghan EOS of order 3.

Good fitting depends on which initial guess parameters
we give.
'''

from math import pi
from numpy import genfromtxt, cos, array
from scipy.optimize import curve_fit
#import matplotlib.pyplot as plt
#from matplotlib import rcParams

##rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True

data = genfromtxt('v.dat', delimiter=' ')
p,v,vca,vc,vo1,vo2 = zip(*data)
p = array(p)
v = array(v)
vca = array(vca)
vc = array(vc)
vo1 = array(vo1)
vo2 = array(vo2)
#vo3 = array(vo3)

def pv_bm3(v,v0,b0,b0p):
    '''
    ======== BM3 p(V) function ==================

    -------- begin parameters -------------------

    v  : volume variable
    v0 : reference volume (equilibrium)
    b0 : static bulk modulus B_0
    b0p: static bulk modulus pressure derivative
         B_0'
    -------- end parameters ---------------------
    '''
    x = v0/v
    f = 0.5*(x**(2./3.) - 1.)
    p = 1.5*b0*f*(2.*f+1.)**(5./2.)*(2.+(3.*(b0p-4.)*f)) 
    return p

def pv_bm3_2(v,v0,b0):
    '''
    ======== BM3 p(V) function ==================

    -------- begin parameters -------------------

    v  : volume variable
    v0 : reference volume (equilibrium)
    b0 : static bulk modulus B_0
    b0p: static bulk modulus pressure derivative
         B_0'
    -------- end parameters ---------------------
    '''
    x = v0/v
    f = 0.5*(x**(2./3.) - 1.)
    p = 3.*b0*f*(2.*f+1.)**(5./2.) 
    return p


# Initial guess vector [a,b,c,d]
p0 = [0.201, 2600., 8.] 
#p0 = [130., 110.] 
#p0 = [1., 1.] 
print vc

# Levenberg-Marquardt
#popt, pcov = curve_fit(pv_bm3, vc, p, p0)
#print 'parameters:', popt

# Evaluation with initial and optimum parameters
f = pv_bm3(vc, *p0)
#f = pv_bm3(p, *popt)
print f
print p

outf = open('fit.out', 'w')
for i in range(len(p)):
    outf.write('{0} {1} {2}\n'.format(vc[i], p[i], f[i]))

#fig = plt.figure()
#ax = fig.add_subplot(1,1,1)

#plt.plot(vca, p, 'o', label='calc')
#plt.plot(vca, f, 
#        label=ur'LJ6-12. \u03B5 ={0:.3f}; \u03C3 ={1:.2f}'.format(*popt))

## Do a legible plot
#ax.set_ylabel(r'E[kcal/mol]', fontsize=20)
#ax.set_xlabel(ur'r[\u00c5]', fontsize=20)
#ax.grid(True)
#plt.yticks(fontsize=14)
#plt.xticks(fontsize=14)
#plt.legend(loc=1)

#plt.savefig('ca.pdf')
