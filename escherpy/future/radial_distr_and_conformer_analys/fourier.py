#!/usr/bin/env python

'''
Trigonometric interpolation with non-linear fitting
Levenberg-Marquardt.

When equally spaced points it is a Discrete Fourier
Transform, but we only want a discrete cosine transform
of n=3 with a particular fitting formula.

Good fitting depends on which initial guess parameters
we give.
'''

from math import pi
from numpy import genfromtxt, cos, array
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


data = genfromtxt('fch2oh.gpl', delimiter=' ')
x,y = zip(*data)
x = array(x)
# degrees -> radians
x = x*pi/180.
# list -> numpy.array
y = array(y)

def fone(w, a, b, n):
    return a*0.5 + 0.5*(b*(1. + cos(n*w)))

def fun(w, a, b, c, d):
    '''
    Discrete cosine transform with n=3
    a) I~naki lecture slides -> parameters represent force constants
    b) (19.13) Arfken 7ed.
    '''
    return a*0.5 + 0.5*(b*(1. + cos(w)) + c*(1. + cos(2*w)) + d*(1. + cos(3*w)))
    #return a*0.5 + 0.5*(b*cos(w) + c*cos(2*w) + d*cos(3*w))


# Initial guess vector [a,b,c,d]
p0 = [-213.81*2., -0.003, -0.003, 0.001]

# Levenberg-Marquardt
popt, pcov = curve_fit(fun, x, y, p0)
print 'parameters', popt

# Evaluation with optimum parameters
f1 = fone(x, popt[0], popt[1], 1.)
f2 = fone(x, popt[0], popt[2], 2.)
f3 = fone(x, popt[0], popt[3], 3.)
f = fun(x, *popt)

# radians -> degrees
x = x*180./pi
plt.plot(x,f1)
plt.plot(x,f2)
plt.plot(x,f3)
plt.plot(x, y, 'o', x, f)
plt.show()






