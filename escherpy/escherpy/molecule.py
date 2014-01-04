# ESCHERpy -- A computational chemistry workflow tool
# Copyright (C) 2013 Daniel Menendez
# 
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

#from __future__ import print_function, division
from math import sqrt, sin, cos, radians, pi
import shlex as lex
from logging import getLogger
import numpy as np
from chem import elements 
from representation import Representation
from grid import Grid
from matplotlib.cm import gray
import matplotlib.pyplot as plt

from atom import colors
from atom import rcov

log = getLogger('escherlog')


class Molecule(object):

    '''
    molecule data structure

    >>> from escherpy import Molecule, escher_data
    >>> mol = Molecule()
    >>> mol.structfile = escher_data + 'mol/ethylene_iso/c2h4.cube'
    >>> mol.readstruct()
    >>> mol.atname
    ['H', 'C', 'H', 'C', 'H', 'H']

    
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

        #: Number of atoms
        self.nat = 0
        ''' number of atoms'''
        #__doc_nat_ = 'number of atoms'
        #: Name of the molecule
        self.name = ""
        ''' name of the molecule'''
        #: Name of the molecule atoms
        self.atname = []
        ''' name of the molecule atoms'''
        #: Atomic numbers
        self.atnumber = np.array([], dtype=np.int)
        ''' atomic numbers'''
        #: Atomic masses
        self.atmass = np.array([])
        ''' atomic masses'''
        #: Atomic coordinates
        self.atxyz = []
        ''' atomic coordinates'''
        self.rep = Representation()
        ''' a :py:class:`Representation` instance'''
        self.rep.window()

        self.cptype = []
        self.cpxyz = np.zeros(3)
        self.vxyz = []


    def openf(self, filename):
        """
        Opens a file and returns lines without whitespaces.

        :param str filename: name of the file 
        :return: a list of the file lines
        :rtype: list of lists

        """
        with open(filename, 'rb') as file:
            text = file.read()

        lines = (line for line in text.splitlines())
        lines = map(lambda s: s.strip(), lines)
        return lines

    def readstruct(self):
        '''
        Reads atomic positions
        '''

        if self.structfile.endswith('.cube'):
            self.readcube()
        if self.structfile.endswith('.xyz'):
            self.readxyz()

    def readcps(self):
        '''
        Reads critical points CPs
        '''

        if self.cpsfile.endswith('CPprop.txt'):
            self.readmultiwfn()

    def readsurf(self, atom):
        '''
        Reads surface vertices
        '''

        if self.basinfile.endswith('.basin'):
            self.readbasin(atom)

    def readbasin(self, atname='H'):
        '''
        Extracts interatomic surface (IAS)  
        from a .basin file
        '''
        log.debug('Reading file {}'.format(self.basinfile))
        lines = self.openf(self.basinfile)

        for n,i in enumerate(lines):

            if n==11:
                # number of vertices
                nv = int(lex.split(i)[0])
                log.debug('Number of vertices being read: {}'.format(nv))
                self.vxyz = np.zeros([nv,3])

            elif n==14:
                for j in xrange(nv):
                    self.vxyz[j,:] = np.array( \
                           [float(k) for k in lex.split(lines[n+j])[0:3]])
                log.debug('Surface vertices coordinates: \n {}'.format(self.vxyz))

        color = [ col/255. for col in colors[atname] ]
        self.rep.surface(nv,self.vxyz, color)
        #for i in xrange(nv):
        #    sphereatom = self.rep.ball(0.05, self.vxyz[i], [1.,1.,1.])
        #    self.rep.ren.AddActor(sphereatom)


    def readmultiwfn(self):
        '''
        Extracts CPs information 
        from a multiwfn critical points file
        '''

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
        '''
        Extracts structure information 
        from a .xyz file
        '''

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
                for j in xrange(self.nat):
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
        '''
        Extracts structure information and grid values
        from a Gaussian .cube file
        '''

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
                for j in xrange(self.nat):
                    self.atnumber = np.append(self.atnumber, 
                            [int(lex.split(lines[n+j])[0])])
                    self.atname.append(elements.find(self.atnumber[-1]))
                    self.atmass = np.append(self.atmass, 
                            [float(lex.split(lines[n+j])[1])])
                    self.atxyz[j,:] = np.array( \
                            [float(k) for k in lex.split(lines[n+j])[2:5]])
                log.debug('Atoms: \n {}'.format(self.atname))
                log.debug('Atomic coordinates: \n {}'.format(self.atxyz))

    def rot(self, *angles):
        '''
        Rotation in 3D 

        :param list angles: [alpha, beta, theta] in degrees
        :returns: R_z(theta) * R_x(alpha) * R_y(beta) * xyz

        '''

        log.debug('Rotating around z axis: {} degrees.'.format(angles[2]))
        log.debug('Rotating around x axis: {} degrees.'.format(angles[0]))
        log.debug('Rotating around y axis: {} degrees.'.format(angles[1]))

        c = [cos(radians(angle)) for angle in angles]
        s = [sin(radians(angle)) for angle in angles]
        rx = np.array([[1,0,0],[0, c[0], -s[0]], [0,s[0],c[0]]])
        ry = np.array([[c[1],0,s[1]],[0, 1, 0], [-s[1],0,c[1]]])
        rz = np.array([[c[2],-s[2],0],[s[2], c[2], 0], [0,0,1]])

        for atom in xrange(len(self.atname)):
            # First: rotate around z axis
            self.atxyz[atom] = np.transpose(np.dot(rz,np.transpose(self.atxyz[atom])))
            # Second: rotate around x axis
            self.atxyz[atom] = np.transpose(np.dot(rx,np.transpose(self.atxyz[atom])))
            # Third: rotate around y axis
            self.atxyz[atom] = np.transpose(np.dot(ry,np.transpose(self.atxyz[atom])))

    def rotx(self, angle):
        '''
        Rotates around the x axis
        '''

        c = cos(radians(angle))
        s = sin(radians(angle))
        rx = np.array([[1,0,0],[0, c, -s], [0,s,c]])

        for atom in xrange(len(self.atname)):
            self.atxyz[atom] = np.transpose(np.dot(rx,np.transpose(self.atxyz[atom])))

    def roty(self, angle):
        '''
        Rotates around the y axis
        '''

        c = cos(radians(angle))
        s = sin(radians(angle))
        ry = np.array([[c,0,s],[0, 1, 0], [-s,0,c]])

        for atom in xrange(len(self.atname)):
            self.atxyz[atom] = np.transpose(np.dot(ry,np.transpose(self.atxyz[atom])))

    def rotz(self, angle):
        '''
        Rotates around the z axis
        '''

        c = cos(radians(angle))
        s = sin(radians(angle))
        rz = np.array([[c,-s,0],[s, c, 0], [0,0,1]])

        for atom in xrange(len(self.atname)):
            self.atxyz[atom] = np.transpose(np.dot(rz,np.transpose(self.atxyz[atom])))


    def stickball(self):
        '''
        Add a stick/ball representation of the structure
        '''
        log.debug('Add stick/ball representation')

        for i in xrange(self.nat):
            color = [ col/255. for col in colors[self.atname[i]] ]
            radi = rcov[self.atname[i]]/200.
            #sphereatom, spheretext = self.rep.ball(radi, self.atxyz[i], color,self.atname[i])
            sphereatom = self.rep.ball(radi, self.atxyz[i], color)
            self.rep.ren.AddActor(sphereatom)
            #self.rep.ren.AddActor(spheretext)
            #spheretext.SetCamera(self.rep.ren.GetActiveCamera())

        for i in xrange(self.nat):
            for j in xrange(i):
                dist = sqrt((self.atxyz[i][0] - self.atxyz[j][0])**2 + 
                            (self.atxyz[i][1] - self.atxyz[j][1])**2 + 
                            (self.atxyz[i][2] - self.atxyz[j][2])**2)
                if (dist) < ((rcov[self.atname[i]]+rcov[self.atname[j]])/50.) :
                    stickbond = self.rep.stick(self.atxyz[i], self.atxyz[j])
                    self.rep.ren.AddActor(stickbond)

    def cpball(self):
        '''
        Add a representation of critical points
        '''
        log.debug('Add critical points ball representation')

        for i in xrange(len(self.cpxyz)):
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

    def isosurface(self, filename, isovalue):
        '''
        Plots an isosurface 

        :param str filename: name of the cube file to read
        :param float isovalue: value of the isosurface to plot

        '''

        self.isovalue = isovalue
        grid = Grid()

        grid.readcube(filename)
        delta = [grid.dx[0,0], grid.dx[1,1], grid.dx[2,2]]

        self.rep.imagedata(grid.x0, delta, grid.n)
        self.rep.outline()
        gridarray = self.rep.array(grid.n, grid.f, 'scalars')
        self.rep.grid.GetPointData().SetScalars(gridarray)
        self.rep.marchingcubes(self.isovalue)

    def nciplot(self, densfile, gradfile):
        '''
        Plots non covalent interactions (NCI) 3D regions
        from dens and grad files generated by critic2

        :param str densfile: Name of the cube file containing 
                             a grid of :math:`\\rho(r)sign(\\lambda_2)`
        :param str gradfile: Name of the cube file containing a grid of the reduced density gradient

        '''

        self.isovalue = 0.5
        rangecolor=[-3.,3.]
        colortable=[[0.,0.,1.],[0.,1.,0.],[1.,0.,0.]]

        self.surfmap(densfile, gradfile, self.isovalue, rangecolor, colortable)

    def nciplot2D(self, densfile, gradfile):
        '''
        Plots non covalent interactions (NCI) in 2D
        from dens and grad files generated by critic2

        :param str densfile: Name of the cube file containing 
                             a grid of :math:`\\rho(r)sign(\\lambda_2)`
        :param str gradfile: Name of the cube file containing a grid of the reduced density gradient

        '''

        dens = Grid()
        dens.readcube(densfile)
        grad = Grid()
        grad.readcube(gradfile)
        #C_s = 1./(2.*(3.*pi**2)**(1./3.))
        #absgrad = np.absolute(grad.alldata)
        #rho4 = np.power(dens.alldata, 4)
        #rho43 = np.power(dens.alldata, -3.)
        #nums = np.multiply(C_s, absgrad)
        #s = np.divide(nums, rho43)
        plt.ylim([-0.,2.])
        plt.xlim([-0.1,0.1])
        dens.alldata = np.divide(dens.alldata, 100.)
        plt.plot(dens.alldata, grad.alldata, 'b8', markersize=3.)
        #tmp = dens.alldata[np.where(dens.alldata < 0.2)]
        #print tmp[np.where(tmp > -0.2)]
        #print dens.alldata[12], grad.alldata[12]
        #plt.scatter(dens.alldata, grad.alldata, c=[1.,0.,0.], s=4.)

        
    def surfmap(self, densfile, gradfile, isovalue=0.5, 
                rangecolor=[-3.,3.],
                colortable=[[0.,0.,1.],[0.,1.,0.],[1.,0.,0.]]):
        '''
        Maps the values of a field at points belonging to a
        surface.
        '''

        self.densfile = densfile
        self.gradfile = gradfile
        self.isovalue = isovalue
        self.rep.lookuptable(rangecolor, colortable)
        self.isosurface(self.gradfile, self.isovalue)

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
        

    def show(self):
        '''
        Renders the molecule
        '''

        self.rep.start()
        plt.show()

                
if __name__ == '__main__':
    import doctest
    doctest.testmod()
                    


