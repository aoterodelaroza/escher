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
        # Jmol colors
        colors = { \
         'H'  : (255,255,255)  ,
         'He' : (217,255,255)  ,
         'Li' : (204,128,255)  ,
         'Be' : (194,255,0)    ,
         'B'  : (255,181,181)  ,
         'C'  : (144,144,144)  ,
         'N'  : (48,80,248)    ,
         'O'  : (255,13,13)    ,
         'F'  : (144,224,80)   ,
         'Ne' : (179,227,245)  ,
         'Na' : (171,92,242)   ,
         'Mg' : (138,255,0)    ,
         'Al' : (191,166,166)  ,
         'Si' : (240,200,160)  ,
         'P'  : (255,128,0)    ,
         'S'  : (255,255,48)   ,
         'Cl' : (31,240,31)    ,
         'Ar' : (128,209,227)  ,
         'K'  : (143,64,212)   ,
         'Ca' : (61,255,0)     ,
         'Sc' : (230,230,230)  ,
         'Ti' : (191,194,199)  ,
         'V'  : (166,166,171)  ,
         'Cr' : (138,153,199)  ,
         'Mn' : (156,122,199)  ,
         'Fe' : (224,102,51)   ,
         'Co' : (240,144,160)  ,
         'Ni' : (80,208,80)    ,
         'Cu' : (200,128,51)   ,
         'Zn' : (125,128,176)  ,
         'Ga' : (194,143,143)  ,
         'Ge' : (102,143,143)  ,
         'As' : (189,128,227)  ,
         'Se' : (255,161,0)    ,
         'Br' : (166,41,41)    ,
         'Kr' : (92,184,209)   ,
         'Rb' : (112,46,176)   ,
         'Sr' : (0,255,0)      ,
         'Y'  : (148,255,255)  ,
         'Zr' : (148,224,224)  ,
         'Nb' : (115,194,201)  ,
         'Mo' : (84,181,181)   ,
         'Tc' : (59,158,158)   ,
         'Ru' : (36,143,143)   ,
         'Rh' : (10,125,140)   ,
         'Pd' : (0,105,133)    ,
         'Ag' : (192,192,192)  ,
         'Cd' : (255,217,143)  ,
         'In' : (166,117,115)  ,
         'Sn' : (102,128,128)  ,
         'Sb' : (158,99,181)   ,
         'Te' : (212,122,0)    ,
         'I'  : (148,0,148)    ,
         'Xe' : (66,158,176)   ,
         'Cs' : (87,23,143)    ,
         'Ba' : (0,201,0)      ,
         'La' : (112,212,255)  ,
         'Ce' : (255,255,199)  ,
         'Pr' : (217,255,199)  ,
         'Nd' : (199,255,199)  ,
         'Pm' : (163,255,199)  ,
         'Sm' : (143,255,199)  ,
         'Eu' : (97,255,199)   ,
         'Gd' : (69,255,199)   ,
         'Tb' : (48,255,199)   ,
         'Dy' : (31,255,199)   ,
         'Ho' : (0,255,156)    ,
         'Er' : (0,230,117)    ,
         'Tm' : (0,212,82)     ,
         'Yb' : (0,191,56)     ,
         'Lu' : (0,171,36)     ,
         'Hf' : (77,194,255)   ,
         'Ta' : (77,166,255)   ,
         'W'  : (33,148,214)   ,
         'Re' : (38,125,171)   ,
         'Os' : (38,102,150)   ,
         'Ir' : (23,84,135)    ,
         'Pt' : (208,208,224)  ,
         'Au' : (255,209,35)   ,
         'Hg' : (184,184,208)  ,
         'Tl' : (166,84,77)    ,
         'Pb' : (87,89,97)     ,
         'Bi' : (158,79,181)   ,
         'Po' : (171,92,0)     ,
         'At' : (117,79,69)    ,
         'Rn' : (66,130,150)   ,
         'Fr' : (66,0,102)     ,
         'Ra' : (0,125,0)      ,
         'Ac' : (112,171,250)  ,
         'Th' : (0,186,255)    ,
         'Pa' : (0,161,255)    ,
         'U'  : (0,143,255)    ,
         'Np' : (0,128,255)    ,
         'Pu' : (0,107,255)    ,
         'Am' : (84,92,242)    ,
         'Cm' : (120,92,227)   ,
         'Bk' : (138,79,227)   ,
         'Cf' : (161,54,212)   ,
         'Es' : (179,31,212)   ,
         'Fm' : (179,31,186)   ,
         'Md' : (179,13,166)   ,
         'No' : (189,13,135)   ,
         'Lr' : (199,0,102)    ,
         'Rf' : (204,0,89)     ,
         'Db' : (209,0,79)     ,
         'Sg' : (217,0,69)     ,
         'Bh' : (224,0,56)     ,
         'Hs' : (230,0,46)     ,
         'Mt' : (235,0,38)     \
        }
        rcov = { \
        'H'  : 37 ,        'Nb' :    137,   'Tl' : 148,
        'He' : 32 ,        'Mo' :    145,   'Pb' : 147,
        'Li' : 134,        'Tc' :    156,   'Bi' : 146,
        'Be' : 90 ,        'Ru' :    126,   'Po' : 200,
        'B'  : 82 ,        'Rh' :    135,   'At' : 200,
        'C'  : 77 ,        'Pd' :    131,   'Rn' : 145,
        'N'  : 75 ,        'Ag' :    153,   'Fr' : 200,
        'O'  : 73 ,        'Cd' :    148,   'Ra' : 200,
        'F'  : 71 ,        'In' :    144,   'Ac' : 200,
        'Ne' : 69 ,        'Sn' :    141,   'Th' : 200,
        'Na' : 154,        'Sb' :    138,   'Pa' : 200,
        'Mg' : 130,        'Te' :    135,   'U'  : 200,
        'Al' : 118,        'I'  :    133,   'Np' : 200,
        'Si' : 111,        'Xe' :    130,   'Pu' : 200,
        'P'  : 106,        'Cs' :    225,   'Am' : 200,
        'S'  : 102,        'Ba' :    198,   'Cm' : 200,
        'Cl' : 99 ,        'La' :    169,   'Bk' : 200,
        'Ar' : 97 ,        'Ce' :    200,   'Cf' : 200,
        'K'  : 196,        'Pr' :    200,   'Es' : 200,
        'Ca' : 174,        'Nd' :    200,   'Fm' : 200,
        'Sc' : 144,        'Pm' :    200,   'Md' : 200,
        'Ti' : 136,        'Sm' :    200,   'No' : 200,
        'V'  : 125,        'Eu' :    200,   'Lr' : 200,
        'Cr' : 127,        'Gd' :    200,   'Rf' : 200,
        'Mn' : 139,        'Tb' :    200,   'Db' : 200,
        'Fe' : 125,        'Dy' :    200,   'Sg' : 200,
        'Co' : 126,        'Ho' :    200,   'Bh' : 200,
        'Ni' : 121,        'Er' :    200,   'Hs' : 200,
        'Cu' : 138,        'Tm' :    200,   'Mt' : 200,
        'Zn' : 131,        'Yb' :    200,          
        'Ga' : 126,        'Lu' :    160,          
        'Ge' : 122,        'Hf' :    150,          
        'As' : 119,        'Ta' :    138,          
        'Se' : 116,        'W'  :    146,          
        'Br' : 114,        'Re' :    159,          
        'Kr' : 110,        'Os' :    128,          
        'Rb' : 211,        'Ir' :    137,          
        'Sr' : 192,        'Pt' :    128,          
        'Y'  : 162,        'Au' :    144,
        'Zr' : 148,        'Hg' :    149,
        }

        for i in range(self.nat):
            # instead use a dictionary
            #if self.atname[i] == 'H':
            #    color = [1,1,1]
            #elif self.atname[i] == 'C':
            #    color = [0.5,0.5,0.5]
            #elif self.atname[i] == 'N':
            #    color = [0,0,1,255]
            #elif self.atname[i] == 'O':
            #    color = [1,0,0]
            #else:
            #    color = [1.0,0.5,0.0]
            color = [ col/255. for col in colors[self.atname[i]]]
            radi = rcov[self.atname[i]]/200.
            sphereatom = self.rep.ball(radi, self.atxyz[i], color)
            self.rep.ren.AddActor(sphereatom)

        for i in range(self.nat):
            for j in range(i):
                dist = sqrt((self.atxyz[i][0] - self.atxyz[j][0])**2 + 
                            (self.atxyz[i][1] - self.atxyz[j][1])**2 + 
                            (self.atxyz[i][2] - self.atxyz[j][2])**2)
                if (dist) < ((rcov[self.atname[i]]+rcov[self.atname[j]])/60.) :
                    stickbond = self.rep.stick(self.atxyz[i], 
                                               self.atxyz[j])
                    self.rep.ren.AddActor(stickbond)

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
    pass
    #import doctest
    #doctest.testmod()
                    


