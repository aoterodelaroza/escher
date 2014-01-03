#!/usr/bin/env python

from math import exp
from random import sample, random
from collections import Counter


T = 300. # K
#k_B = 1.3806488e-23
R = 8.3144621 # J K^{-1}mol^{-1}

energies = [0.0, 0.2, 1.5, 2.9, 7.3] # KJ mol^{-1}
for i in range(len(energies)):
	energies[i] = energies[i]*1000. # J mol^{-1}

dipoles = [0.0, 1.2, 1.1, 1.7, 2.9] # D

beta = 1./(T*R)

print '======================'
print 'Boltzmann distribution'
print '======================'
print ' State  Probability'
print '-------------------'
Z = 0.
for i in range(5):
	Z = Z + exp(-beta*energies[i])

P = []
for i in range(5):
	P.append(exp(-beta*energies[i])/Z)
	print i, P[i]

print '======================'
print 'Monte-Carlo trajectory'
print '======================'
v = []
# Choose amongst the posible states
v.append(sample(range(1,6),1)[0] - 1)
for n in range(1000000):
	v.append(sample(range(1,6),1)[0] - 1)
	ae = energies[v[n+1]] - energies[v[n]]
	if ae>0.0:
		if exp(-beta*ae)<random():
			v[n+1] = v[n]

# Count occurrences of the same value
c = Counter(v)

p = []
dipole = 0.0
print ' State  Population  Probability'
print '-------------------------------'
for i in range(5):
	p.append(float(c[i])/float(len(v)))
	print i, c[i], p[i]
	dipole = dipole + p[i]*dipoles[i]
	
dipole = dipole/5.
print
print 'Dipole: ', dipole, 'D'
	






