#!/usr/bin/env python

from numpy import linspace,exp
from numpy import genfromtxt
from numpy.random import randn
from scipy.interpolate import UnivariateSpline
from scipy.integrate import quad
import matplotlib.pyplot as plt


data = genfromtxt('rdf_na_wat.csv', delimiter=',')
x,y = zip(*data)

s = UnivariateSpline(x, y, s=0)
xs = linspace(2.13, 3.1, 1000000)
ys = s(xs)

print('Integral of g(r) is:')
print(s.integral(2.0,3.1))

def fun(x):
    return s(x)*(x**2)


res = quad(fun, 2.13,3.1)
print('Integral of g(r)r^2 is:')
print('Result              Error')
print(res)


plt.plot(x,y)
plt.plot(xs,ys)
plt.show()






# = linspace(-3, 3, 100)
# = exp(-x**2) + randn(100)/10
#
# = UnivariateSpline(x, y, s=1)
#s = linspace(-3, 3, 1000)
#ys = s(xs)
