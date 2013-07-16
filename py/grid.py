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

from logging import getLogger
import shlex as lex
import numpy as np
log = getLogger('escherlog')

class Grid():

    '''
    function g = grid()

    grid - create an empty grid structure.

    Output:
    {g}: the empty grid structure with all the fields defined.
    '''

    def __init__(self):

        self.name = ""
        self.x0 = np.matrix([0,0,0])
        self.dx = np.zeros([3,3])
        self.a = np.zeros([3,3])
        # cube grid dimensions
        self.n = np.array([], dtype=np.int)
        # all cube scalr data dim=1 array
        self.alldata = []
        self.f = []
        self.omega = 0.0

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
                self.x0 = [float(j) for j in lex.split(i)[1:4]]
                log.debug('Grid origin {0:.6f} {1:.6f} {2:.6f}'.format(*self.x0))
            elif n==3:
                for j in range(3):
                    # reading grid dimensions
                    self.n = np.append(self.n, [int(lex.split(lines[n+j])[0])])
                    # reading grid spacings
                    self.dx[j,:] = np.array([float(k) for k in lex.split(lines[n+j])[1:4]])
                    log.debug('Grid matrix {0:15.9f} {1:15.9f} {2:15.9f}'.format(*self.dx[j,:]))
                    self.a[j,:] = self.dx[j,:].dot(self.n[j])
                log.debug('Grid dimensions {0:d} {1:d} {2:d}'.format(*self.n))
                self.omega = np.linalg.det(self.a)
            elif n==(6+self.nat):
                for j in range(len(lines)-6-self.nat):
                    self.alldata.extend([float(k) for k in lex.split(lines[n+j])])
                self.alldata = np.array(self.alldata)
                log.debug('Grid data points {0:d}'.format(len(self.alldata)))
        self.f = np.zeros(self.n)
        ndim = 0
        for i in range(self.n[0]):
            for j in range(self.n[1]):
                self.f[i,j,:] = self.alldata[ndim:ndim+self.n[2]] 
                ndim = ndim + self.n[2]


