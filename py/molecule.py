# Copyright (C) 2011 Victor Lua~na and Alberto Otero-de-la-Roza
#
# This octave routine is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version. See <http://www.gnu.org/licenses/>.
#
# The routine distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.

#from __future__ import print_function, division
#from math import cos, sin, pi, sqrt, acos, degrees
from math import sqrt
import shlex as lex
import numpy as np
#import pybel
from chem import elements 
import logcolors
from representation import Representation
from grid import Grid

log = logcolors.escherlog()
#log.setLevel(log.DEBUG)

class Molecule():

    '''
    molecule data structure
    >>> mol = Molecule()
    >>> mol.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
    >>> mol.a
    array([[ 13.33332 ,   0.      ,   0.      ],
           [  0.      ,  13.33332 ,   0.      ],
           [  0.      ,   0.      ,  15.999984]])

    
    '''

    def __init__(self):
        '''
        molecule - create an empty molecule structure and initialize Molecule().nat
        '''
        #log.critical('critical message example')
        #log.debug('debug message example')
        #log.warning('warning message example')
        #log.error('error message example')

        log.debug('Initializing data values')

        # Number of atoms
        self.nat = 0
        ''' number of atoms'''
        __doc_nat_ = 'number of atoms'
        # Name of the molecule
        self.name = ""
        # Name of the molecule atoms
        self.atname = []
        # Atomic numbers
        self.atnumber = np.array([], dtype=np.int)
        # Atomic masses
        self.atmass = np.array([])
        # Atomic coordinates
        self.atxyz = []
#        # metric matrix
#        self.a = np.zeros([3,3])
#        # cube grid spacings
#        self.dx = self.a
#        # cube grid dimensions
#        self.n = np.array([], dtype=np.int)
#        # all cube scalr data dim=1 array
#        self.alldata = []


    def abrir(self, filename):
        """
        Opens a file and returns lines without whitespaces.
        """
        with open(filename, 'r') as file:
            text = file.read()

        lines = [line for line in text.splitlines() if line]
        lines = map(lambda s: s.strip(), lines)
        return lines

    def readcube(self, filename):
        """
        Extracts structure information and grid values
        from a Gaussian *.cube file
        """

        ## Openbabel data 
        #mol = pybel.readfile("cube", filename).next()
        #file = open(filename, 'r')
        #self.name = file.readline()
        #self.title = file.readline()
        #self.nat = file.readline()[0]
        #self.nat = len(mol.atoms)
        #for atom in mol:
            #self.atnumber.append(atom.atomicnum)
            #self.atname.append(elements.find(atom.atomicnum))
            #self.atmass.append(atom.atomicmass)
            #self.atxyz.append(atom.coords)
        log.debug('Reading file {}'.format(filename))

        ## parsing is not done the clever way, it could be improved
        lines = self.abrir(filename)
        for n,i in enumerate(lines):
            if n==0: # first line
                self.name = i
            elif n==1: # second line
                self.title = i
            elif n==2: # third line
                self.nat = int(lex.split(i)[0])
#                self.x0 = [float(j) for j in lex.split(i)[1:4]]
#                log.info('Grid origin {0:.6f} {1:.6f} {2:.6f}'.format(*self.x0))
                # initialize atomic names, coordinates with right size
                self.atxyz = np.zeros([self.nat,3])
#            elif n==3:
#                for j in range(3):
#                    # reading grid dimensions
#                    self.n = np.append(self.n, [int(lex.split(lines[n+j])[0])])
#                    # reading grid spacings
#                    self.dx[j,:] = np.array([float(k) for k in lex.split(lines[n+j])[1:4]])
#                    log.info('Grid matrix {0:15.9f} {1:15.9f} {2:15.9f}'.format(*self.dx[j,:]))
#                    self.a[j,:] = self.dx[j,:].dot(self.n[j])
#                log.info('Grid dimensions {0:d} {1:d} {2:d}'.format(*self.n))
#                self.omega = np.linalg.det(self.a)
            elif n==6:
                for j in range(self.nat):
                    self.atnumber = np.append(self.atnumber, [int(lex.split(lines[n+j])[0])])
                    self.atname.append(elements.find(self.atnumber[-1]))
                    self.atmass = np.append(self.atmass, [float(lex.split(lines[n+j])[1])])
                    self.atxyz[j,:] = np.array([float(k) for k in lex.split(lines[n+j])[2:5]])
#            elif n==(6+self.nat):
#                for j in range(len(lines)-6-self.nat):
#                    self.alldata.extend([float(k) for k in lex.split(lines[n+j])])
#                self.alldata = np.array(self.alldata)
#                log.info('Grid data points {0:d}'.format(len(self.alldata)))
#        self.f = np.zeros(self.n)
#        ndim = 0
#        for i in range(self.n[0]):
#            for j in range(self.n[1]):
#                self.f[i,j,:] = self.alldata[ndim:ndim+self.n[2]] 
#                ndim = ndim + self.n[2]

    def nciplot(self, file1, file2):
        densgrid = Grid()
        gradgrid = Grid()

        #densgrid.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')

    def stickball(self):
        rep = Representation()
        rep.window()
        log.debug('Add stick/ball representation')
        for i in range(self.nat):
            # instead use a dictionary
            if self.atname[i] == 'H':
                color = [1,1,1]
            if self.atname[i] == 'C':
                color = [0.5,0.5,0.5]
            if self.atname[i] == 'N':
                color = [0,0,1,255]
            if self.atname[i] == 'O':
                color = [1,0,0]
            else:
                color = [0.5,0.5,0.5]
            sphereatom = rep.ball(0.6, self.atxyz[i], color)
            rep.ren.AddActor(sphereatom)

        for i in range(self.nat):
            for j in range(i):
                dist = sqrt((self.atxyz[i][0] - self.atxyz[j][0])**2 + 
                            (self.atxyz[i][1] - self.atxyz[j][1])**2 + 
                            (self.atxyz[i][2] - self.atxyz[j][2])**2)
                if (dist) < 3. :
                    rep.ren.AddActor(rep.stick(self.atxyz[i], 
                                      self.atxyz[j]))
        self.rep = rep

    def isosurface(self):
        dens = Grid()
        dens.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
        #delta = [dens.dx[0,0]/dens.n[0], dens.dx[1,1]/dens.n[1], 
        #         dens.dx[2,2]/dens.n[2]]
        delta = [dens.dx[0,0], dens.dx[1,1], dens.dx[2,2]]
        self.rep.gridvolume(dens.x0, delta, dens.n, dens.f)
        self.rep.marchingcubes()
        
    def start(self):

        self.rep.start()


                
if __name__ == '__main__':
    import doctest
    doctest.testmod()
                    


