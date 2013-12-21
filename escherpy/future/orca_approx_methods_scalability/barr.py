#!/usr/bin/env python

'''
Data analysis of RI, RIJCOSX
with MP2 within water clusters
'''

from __future__ import print_function, division
import os
import shlex as sh
import matplotlib.pyplot as plt
#from matplotlib import rcParams
import numpy as np
from scipy.optimize import curve_fit

#rcParams['text.usetex']=True
#rcParams['text.latex.unicode']=True

def f(m, a, b):
    '''
    Function used to fit data
    '''
    return a*np.power(m,b)

def extract(calc, clusters):
    '''
    Extracts data for a given type of calculation
    '''
    os.remove('{0}.txt'.format(calc))
    print('M', 'time', file=open('{0}.txt'.format(calc), 'a'))
    nbas = []
    t = []
    energy = []
    for i in clusters:
        inf = open('{0}_{1}.out'.format(i, calc), 'r')
        lines = inf.readlines()
        for j in lines:
            if 'TOTAL RUN TIME' in j:
                last = sh.split(j)
                time = float(last[3])*24.*3600. \
                     + float(last[5])*3600.     \
                     + float(last[7])*60.       \
                     + float(last[9])           \
                     + float(last[11])*0.001
            #if 'Dimension of the basis' in j:
            #if 'of primitive gaussian functions' in j:
            if 'of contracted basis functions' in j:
                nbasis = int(sh.split(j)[-1])
            if 'FINAL SINGLE POINT ENERGY' in j:
                e = float(sh.split(j)[-1])
        print(nbasis, time, file=open('{0}.txt'.format(calc), 'a'))
        nbas.append(nbasis)
        t.append(time)
        energy.append(e)
    return nbas, t, energy

def fit(calc):
    '''
    Fits and plots data for a given type of calculation
    '''

    nbas, t, e = extract(calc, clusters)
    energy.append(e)

    x = np.linspace(nbas[0], nbas[-1], 100)
    nbas = np.array(nbas)
    t = np.array(t)

    # Initial guess vector of parameters
    p0 = [0.0001, 2.59]
    # Levenberg-Marquardt non linear fit
    popt, pcov = curve_fit(f, nbas, t, p0)
    # Values calculated with the optimum parameters
    pf = f(x, *popt)

    print('=========')
    print(calc)
    print('=========')
    print('prefactor', popt[0])
    print('exponent',  popt[1])
    plt.plot(nbas, t, 'o', label='{0} data'.format(calc))
    plt.plot(x, pf, 
        label='{0} fit. y = {1:.2e}M^{2:5.2f}'.format(calc, *popt))
    return energy

clusters = [2, 3, 8, 9]
energy = []

# Create plot window
fig = plt.figure()
ax = fig.add_subplot(1,1,1)

energy = fit('mp2')
energy = fit('rimp2')
energy = fit('rijcosx')
energy = np.array(energy)
energy = energy.T
energy = energy*627.509391

# Errors
print('RI error', 'RIJCOSX error')
for i,j in zip(energy, clusters):
    print(abs(i[0]-i[1])/(j*3.),abs(i[0]-i[2])/(j*3.))

# Do a legible plot
ax.set_ylabel(r'Time[s]', fontsize=26)
ax.set_xlabel(r'M', fontsize=26)
ax.grid(True)
plt.yticks(fontsize=16)
plt.xticks(fontsize=16)
plt.legend(loc=2)

plt.savefig('fit.pdf', dpi=300)

