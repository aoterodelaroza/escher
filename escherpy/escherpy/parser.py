#!/usr/bin/env python

import os
import csv
import subprocess
import numpy as np
import shlex as lex
from logging import getLogger
from chem import elements 
#from representation import Representation

log = getLogger('escherlog')


class Parser(object):


    def __init__(self):

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
        #self.rep = Representation()
        #''' a :py:class:`Representation` instance'''
        #self.rep.window()
        self.n = np.array([], dtype=np.int)
        '''cube grid dimensions'''
        #self.x0 = np.matrix([0,0,0])
        self.dx = np.zeros([3,3])
        self.a = np.zeros([3,3])
        #self.f = []
        #self.omega = 0.0
        self.alldata = []

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

    def critic2(self):
        """
        Extracts basin charges, ...
        from a .outcritic file
        """
        lines = self.openf(filename)

        charges = []
        volumes = []

        for n,i in enumerate(lines):

            if i.startswith('Number of atoms in the unit cell'):
                self.nat = int(i.split()[-1])

            if i.startswith('* List of basins and local properties'):
                for j in range(2,self.nat+1):
                    charges.append(float(lines[j+n].split()[-1]))
                    volumes.append(float(lines[j+n].split()[-2]))
                self.charges = charges
                self.volume = volumes

            if i.startswith('Total integration'):
                self.totcharge = float(i.split()[-1])

    def multiwfn(self,filename):
        '''
        Extracts CPs information 
        from a multiwfn critical points file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)
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


    def qe(self):
        ######################################################
        ### SCF
        ######################################################
        """
        Extracts data from {prefix}.scf.out
        """
        log.info('Extracting from {}.scf.out...'.format(self.prefix))
        lines = self.abrir('{}.scf.out'.format(self.prefix))
        ks = {}
        for n,i in enumerate(lines):
            if 'Quantum ESPRESSO suite' in i:
                self.espresso = True
                log.info('Quantum ESPRESSO calculation.')
            if i.startswith('JOB DONE'):
                self.sucesscf = True
                log.info('JOB DONE calculation.')
            if i.startswith('convergence has been achieved'):
                self.niter = int(sh.split(i)[-2])
                log.info('Converged.')

            if i.startswith('!'):
                self.energies = []
                self.energies.append(float(sh.split(i)[-2]))
                self.energies_units = sh.split(i)[-1]
                log.debug('Energy.')
            if 'bravais-lattice' in i:
                self.ibrav = int(sh.split(i)[-1])
            if i.startswith('PWSCF'):
                times = [sh.split(i)[x] for x in [2, 4]]
                self.times = map(lambda s : float(s.strip('s')), times)
            if i.startswith('unit-cell volume'):
                self.volume = float(sh.split(i)[-2])
                self.volume_units = sh.split(i)[-1]
            if i.startswith('number of atoms/cell'):
                self.natom = int(i.rsplit(' ', 1)[1])
            if i.startswith('number of atomic types'):
                self.ntype = int(i.rsplit(' ', 1)[1])
            if i.startswith('number of electrons'):
                self.nelec = float(i.rsplit(' ', 1)[1])

            # Cell parameters
            if i.startswith('lattice parameter'):
                self.alat = float(sh.split(i)[-2])
                self.alat_units = sh.split(i)[-1]
            if i.startswith('celldm(1)'):
                alats = [float(sh.split(i)[x]) for x in [1,3,5]]
            if i.startswith('celldm(4)'):
                alatsd = [float(sh.split(i)[x]) for x in [1,3,5]]

            # Metrics matrix
            if i.startswith('a(1)'):
                v1 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('a(2)'):
                v2 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('a(3'):
                v3 = [float(sh.split(i)[x]) for x in [3,4,-2]]

            # Metrics inverse matrix
            if i.startswith('b(1)'):
                w1 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('b(2)'):
                w2 = [float(sh.split(i)[x]) for x in [3,4,-2]]
            if i.startswith('b(3'):
                w3 = [float(sh.split(i)[x]) for x in [3,4,-2]]

            # K points
            if i.startswith('number of k points'):
                self.nkpoints = int(sh.split(i)[-1])
                n = n + 1
                for k in range(self.nkpoints):
                    n = n + 1
                    ks['{}'.format(k)] = \
                        [sh.split(lines[n])[y] for y in [-6, -5, -4, -1]]
                    ks['{}'.format(k)] = \
                            map(lambda s : float(s.strip('),')),
                                ks['{}'.format(k)])
                    self.ks = ks
                log.debug('K points.')

        self.celldm = alats + alatsd
        log.debug('Cell parameters.')
        self.R = np.array([v1,v2,v3])
        log.debug('Metrics matrix.')
        self.rR = np.array([v1,v2,v3])
        log.debug('Metrics inverse matrix.')

        ######################################################
        ### DOS
        ######################################################
        """
        Extracts data form {prefix}.dos.out
        """
        dosout = '{}.dos.out'.format(self.prefix)
        log.info('Extracting from: {}...'.format(dosout))
        if not os.path.isfile(dosout):
            return 'No efermi taken still. Run dosbands(), then plotbandx()'

        lines = self.abrir(dosout)
        for i in lines:
            if i.startswith('the Fermi energy is '):
                self.efermi = float(sh.split(i)[-2])
                self.efermi_units = sh.split(i)[-1]
                log.debug('Efermi.')

        pass

    def xyz(self, filename):
        '''
        Extracts structure information 
        from a .xyz file
        '''

        log.debug('Reading file {}'.format(filename))

        ## parsing is not done the clever way, it could be improved
        lines = self.openf(filename)
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

    def cubeheader(self, filename):
        '''
        Extracts structure information and grid values
        from a Gaussian .cube file
        '''

        log.debug('Reading file {}'.format(filename))

        ## parsing is not done the clever way, it could be improved
        lines = self.openf(filename)
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

    def cube(self, filename):
        """
        Extracts structure information and grid values
        from a Gaussian *.cube file
        """

        log.debug('Reading file {}'.format(filename))

        ## parsing is not done the clever way, it could be improved
        lines = self.openf(filename)
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
                #self.alldata = np.genfromtxt(lines[n:])
                for j in csv.reader(lines[n:], delimiter=' '):
                    j = filter(None, j)
                    j = map(float, j)
                    self.alldata.extend(j)
                #for j in range(len(lines)-6-self.nat):
                    #self.alldata.extend([float(k) for k in lex.split(lines[n+j])])
                    #self.alldata = np.append(self.alldata, np.fromstring(lines[n+j], sep=' '))
                self.alldata = np.array(self.alldata)
                log.debug('Grid data points {0:d}'.format(len(self.alldata)))
        self.f = np.zeros(self.n)
        ndim = 0
        for i in range(self.n[0]):
            for j in range(self.n[1]):
                self.f[i,j,:] = self.alldata[ndim:ndim+self.n[2]] 
                ndim = ndim + self.n[2]


    def basin(self, filename):
        '''
        Extracts interatomic surface (IAS)  
        from a .basin file
        '''
        log.debug('Reading file {}'.format(filename))
        lines = self.openf(filename)

        for n,i in enumerate(lines):

            if n==11:
                # number of vertices
                nv = int(lex.split(i)[0])
                log.debug('Number of vertices being read: {}'.format(nv))
                vxyz = np.zeros([nv,3])

            elif n==14:
                for j in xrange(nv):
                    vxyz[j,:] = np.array( 
                           [float(k) for k in lex.split(lines[n+j])[0:3]])
                log.debug('Surface vertices coordinates: \n {}'.format(
                                                            vxyz))
        self.nv = nv
        self.vxyz = vxyz


    def pmd(self):
        '''
        Extracts data from a Promolden run
        .pmdout
        '''
        log.debug('Reading file {}'.format(filename))
        lines = self.openf(filename)

    def henkelman(self, filename):
        '''
        Extracts data from a Henkelman's bader run
        ACF.dat
        '''
        log.debug('Reading file {}'.format(filename))
        lines = self.openf(filename)

    def orca(self):
        '''
        Extracts information 
        from a orca output file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)

    def nwchem(self):
        '''
        Extracts information 
        from a NWChem output file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)

    def g09(self):
        '''
        Extracts information 
        from a g09 output file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)

    def gamess(self):
        '''
        Extracts information 
        from a gamess output file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)

    def vasp(self):
        '''
        Extracts information 
        from a VASP output file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)

    def abinit(self):
        '''
        Extracts information 
        from a abinit output file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)

    def psi4(self):
        '''
        Extracts information 
        from a PSI4 output file
        '''

        log.debug('Reading file {}'.format(filename))

        lines = self.openf(filename)



