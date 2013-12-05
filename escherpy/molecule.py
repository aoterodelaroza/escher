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
from math import sqrt
import shlex as lex
from logging import getLogger
import numpy as np
from chem import elements 
from representation import Representation
from grid import Grid

from atom import colors
from atom import rcov

log = getLogger('escherlog')


#log = logcolors.escherlog()
#log.setLevel(log.DEBUG)

class Molecule():

    '''
    molecule data structure
    >>> mol = Molecule()
    >>> mol.structfile = '../../escher_data/mol/ethylene_iso/c2h4.cube'
    >>> mol.readstruct()
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
        self.rep = Representation()
        self.rep.window()

        self.cptype = []
        self.cpxyz = np.zeros(3)
        self.vxyz = []


    def openf(self, filename):
        """
        Opens a file and returns lines without whitespaces.
        """
        with open(filename, 'r') as file:
            text = file.read()

        lines = [line for line in text.splitlines()]
        lines = map(lambda s: s.strip(), lines)
        return lines

    def readstruct(self):

        if self.structfile.endswith('.cube'):
            self.readcube()
        if self.structfile.endswith('.xyz'):
            self.readxyz()

    def readcps(self):

        if self.cpsfile == 'CPprop.txt':
            self.readmultiwfn()

    def readsurf(self, atom):

        if self.basinfile.endswith('.basin'):
            self.readbasin(atom)

    def readbasin(self, atname='H'):
        """
        Extracts interatomic surface (IAS)  
        from a *.basin file
        """
        log.debug('Reading file {}'.format(self.basinfile))
        lines = self.openf(self.basinfile)

        for n,i in enumerate(lines):

            if n==11:
                # number of vertices
                nv = int(lex.split(i)[0])
                log.debug('Number of vertices being read: {}'.format(nv))
                self.vxyz = np.zeros([nv,3])

            elif n==14:
                for j in range(nv):
                    self.vxyz[j,:] = np.array( \
                           [float(k) for k in lex.split(lines[n+j])[0:3]])
                log.debug('Surface vertices coordinates: \n {}'.format(self.vxyz))

        color = [ col/255. for col in colors[atname] ]
        self.rep.surface(nv,self.vxyz, color)
        #for i in range(nv):
        #    sphereatom = self.rep.ball(0.05, self.vxyz[i], [1.,1.,1.])
        #    self.rep.ren.AddActor(sphereatom)


    def readmultiwfn(self):
        """
        Extracts CPs information 
        from a multiwfn critical points file
        """

        log.debug('Reading file {}'.format(self.cpsfile))

        lines = self.openf(self.cpsfile)
        for n,i in enumerate(lines):

            if i.startswith('======='): 

                cptype = lex.split(i)[4]
                if cptype == '(3,-3)':
                    # nuclear critical point
                    self.cptype = np.append(self.cptype, 'NCP')
                if cptype == '(3,-1)':
                    # bond critical point
                    self.cptype = np.append(self.cptype, 'BCP')
                if cptype == '(3,1)':
                    # ring critical point
                    self.cptype = np.append(self.ctypee, 'RCP')
                if cptype == '(3,3)':
                    # cage critical point
                    self.cptype = np.append(self.cptype, 'CCP')


            if i.startswith('Position'): 

                cpx = float(lex.split(i)[-3])
                cpy = float(lex.split(i)[-2])
                cpz = float(lex.split(i)[-1])
                self.cpxyz = np.vstack([self.cpxyz, [cpx,cpy,cpz]])

        log.debug('Critical points type: \n {}'.format(self.cptype))
        # Delete first row of the matrix, the one used to initialize
        self.cpxyz = np.delete(self.cpxyz, 0, 0)
        log.debug('Critical points coordinates: \n {}'.format(self.cpxyz))

    def readxyz(self):
        """
        Extracts structure information 
        from a *.xyz file
        """

        log.debug('Reading file {}'.format(self.structfile))

        ## parsing is not done the clever way, it could be improved
        lines = self.openf(self.structfile)
        for n,i in enumerate(lines):
            if n==0: # first line
                self.nat = int(lex.split(i)[0])
                log.debug('Number of atoms: {}'.format(self.nat))
                # initialize atomic names, coordinates with right size
                self.atxyz = np.zeros([self.nat,3])
            elif n==1: # second line
                self.title = i
                log.debug('Title: {}'.format(self.title))
            elif n==2:
                for j in range(self.nat):
                    self.atname = np.append(self.atname, 
                            [str(lex.split(lines[n+j])[0])])
                    #self.atname.append(elements.find(self.atnumber[-1]))
                    #self.atmass = np.append(self.atmass, 
                    #        [float(lex.split(lines[n+j])[1])])
                    self.atxyz[j,:] = np.array( \
                            [float(k) for k in lex.split(lines[n+j])[1:4]])
                log.debug('Atoms: \n {}'.format(self.atname))
                log.debug('Atomic coordinates: \n {}'.format(self.atxyz))

    def readcube(self):
        """
        Extracts structure information and grid values
        from a Gaussian *.cube file
        """

        log.debug('Reading file {}'.format(self.structfile))

        ## parsing is not done the clever way, it could be improved
        lines = self.openf(self.structfile)
        for n,i in enumerate(lines):
            if n==0: # first line
                self.name = i
            elif n==1: # second line
                self.title = i
            elif n==2: # third line
                self.nat = int(lex.split(i)[0])
                log.debug('Number of atoms: {}'.format(self.nat))
                # initialize atomic names, coordinates with right size
                self.atxyz = np.zeros([self.nat,3])
            elif n==6:
                for j in range(self.nat):
                    self.atnumber = np.append(self.atnumber, 
                            [int(lex.split(lines[n+j])[0])])
                    self.atname.append(elements.find(self.atnumber[-1]))
                    self.atmass = np.append(self.atmass, 
                            [float(lex.split(lines[n+j])[1])])
                    self.atxyz[j,:] = np.array( \
                            [float(k) for k in lex.split(lines[n+j])[2:5]])
                log.debug('Atoms: \n {}'.format(self.atname))
                log.debug('Atomic coordinates: \n {}'.format(self.atxyz))

    def stickball(self):
        log.debug('Add stick/ball representation')

        for i in range(self.nat):
            color = [ col/255. for col in colors[self.atname[i]] ]
            radi = rcov[self.atname[i]]/200.
            #sphereatom, spheretext = self.rep.ball(radi, self.atxyz[i], color,self.atname[i])
            sphereatom = self.rep.ball(radi, self.atxyz[i], color)
            self.rep.ren.AddActor(sphereatom)
            #self.rep.ren.AddActor(spheretext)
            #spheretext.SetCamera(self.rep.ren.GetActiveCamera())

        for i in range(self.nat):
            for j in range(i):
                dist = sqrt((self.atxyz[i][0] - self.atxyz[j][0])**2 + 
                            (self.atxyz[i][1] - self.atxyz[j][1])**2 + 
                            (self.atxyz[i][2] - self.atxyz[j][2])**2)
                if (dist) < ((rcov[self.atname[i]]+rcov[self.atname[j]])/50.) :
                    stickbond = self.rep.stick(self.atxyz[i], 
                                               self.atxyz[j])
                    self.rep.ren.AddActor(stickbond)

    def cpball(self):
        log.debug('Add critical points ball representation')

        for i in range(len(self.cpxyz)):
            if self.cptype[i] == 'NCP':
                color = [0.2,0.2,0.2]
                radi = 0.0
            if self.cptype[i] == 'BCP':
                color = [0.1,0.95,0.1]
                radi = 0.2
            if self.cptype[i] == 'RCP':
                color = [1,1,1]
                radi = 0.1
            if self.cptype[i] == 'CCP':
                color = [1,1,1]
                radi = 0.1

            #sphereatom, spheretext = self.rep.ball(radi, self.atxyz[i], color,self.atname[i])
            sphereatom = self.rep.ball(radi, self.cpxyz[i], color)
            self.rep.ren.AddActor(sphereatom)
            #self.rep.ren.AddActor(spheretext)
            #spheretext.SetCamera(self.rep.ren.GetActiveCamera())

    def isosurface(self, filename):

        grid = Grid()

        grid.readcube(filename)
        delta = [grid.dx[0,0], grid.dx[1,1], grid.dx[2,2]]

        self.rep.imagedata(grid.x0, delta, grid.n)
        self.rep.outline()
        gridarray = self.rep.array(grid.n, grid.f, 'scalars')
        self.rep.grid.GetPointData().SetScalars(gridarray)
        self.rep.marchingcubes(self.isovalue)
        
    def nciplot(self):

        self.isovalue = 0.5
        self.rep.lookuptable()
        self.isosurface(self.gradfile)

        dens = Grid()
        dens.readcube(self.densfile)
        densarray = self.rep.array(dens.n, dens.f, 'dens')
        self.rep.grid.GetPointData().AddArray(densarray)

        self.rep.mapper.ScalarVisibilityOn()
        self.rep.mapper.SetLookupTable(self.rep.colorNCI)
        self.rep.mapper.SetScalarModeToUsePointFieldData()
        self.rep.mapper.SelectColorArray(1)
        self.rep.mapper.SetScalarRange(-5.0, 5.0)
        self.rep.scalarbar()
        

    def start(self):

        self.rep.start()

                
if __name__ == '__main__':
    import doctest
    doctest.testmod()
                    


