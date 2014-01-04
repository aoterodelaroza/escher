#!/usr/bin/env python

import os
import glob
import pexpect
import subprocess
import numpy as np
from logging import getLogger
from templates.orca import orcainp
from templates.multiwfn import multiwfnin

log = getLogger('escherlog')

class Orca(object):

    def __init__(self):

        self.ang_to_au = 1.0 / 0.5291772083

        self.elements = [None,
             "H", "He",
             "Li", "Be",
             "B", "C", "N", "O", "F", "Ne",
             "Na", "Mg",
             "Al", "Si", "P", "S", "Cl", "Ar",
             "K", "Ca",
             "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
             "Ga", "Ge", "As", "Se", "Br", "Kr",
             "Rb", "Sr",
             "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd",
             "In", "Sn", "Sb", "Te", "I", "Xe",
             "Cs", "Ba",
             "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
             "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg",
             "Tl", "Pb", "Bi", "Po", "At", "Rn",
             "Fr", "Ra",
             "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm", "Md", "No",
             "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Uub"]

    def _run(self, basename):

        orcadict = {'basename': basename}

        orcainpf = open(basename + ".inp", "w")
        orcainpf.write(orcainp.format(**orcadict))

        #orcaout = subprocess.check_output(['orca' , basename + '.inp'])
        #with open(basename + ".out", 'wb') as orcaoutf:
        #    orcaoutf.write(orcaout)

    def run(self, basename):
        '''
        Requires: .xyz 

        Output: .gbw .out .scfp
        '''

        self._run(basename)
        orcaout = subprocess.check_output(['orca' , basename + '.inp'])
        with open(basename + ".out", 'wb') as orcaoutf:
            orcaoutf.write(orcaout)

        os.remove(basename + '.prop')

    def genwfn(self,basename):
        '''
        Generates a .wfn file from a Orca calculation .gbw

        Requires a molden2aim executable properly added to the PATH
        '''

        subprocess.check_call(["orca_2mkl", basename , "-molden"])

        with open(basename + ".molden", 'wb') as m2ain:
            m2ain.write('[Program] Orca\n')
            for line in open(basename + '.molden.input', 'rb'):
                m2ain.write(line)

        child = pexpect.spawn('molden2aim')

        child.logfile = open('molden2aim.log.tmp', 'w')

        child.expect('.*Press <ENTER> to continue\r\n')
        child.sendline('')

        child.expect(' > ')
        child.sendline(basename + '.molden \r\n')
        child.sendline('')

        child.expect(' > ')
        child.sendline('')

        child.expect('.*Press <ENTER> to exit\r\n')
        child.sendline('')
        
        os.remove(basename + '.molden')
        os.remove(basename + '_new.molden')
        os.remove(basename + '.molden.input')
        os.remove('molden2aim.log.tmp')


    def read_xyz(self, xyz):

        atoms = []
        x = [] 
        y = []
        z = []
        f = open(xyz, "r")

        nat = int(f.readline().strip())
        f.next()

        for n,line in enumerate(f):

            if n == nat:
                break

            data = line.split()
            atoms.append(data[0])
            x.append(float(data[1]))
            y.append(float(data[2]))
            z.append(float(data[3]))

        f.close()

        return atoms, np.array(x), np.array(y), np.array(z)


    def read_vpot(self, vpot):

        v = []
        f = open(vpot, "r")
        f.next()

        for line in f:
            data = line.split()
            v.append(float(data[3]))
        f.close()

        return np.array(v)

    def potcube(self, basename, npoints=80, margin=2):
        '''
        Create a MEP .cube file with ORCA.

        @ basename - file name without the extension;
                     this should be the same for the .gbw and .scfp.
        @ npoints  - number of grid points per original cell size
                     (80 should be fine)

        Dependencies: .xyz .gbw .scfp
        orca_vpot has to be properly added to the PATH.
        '''

        npoints = npoints*margin

        atoms, x, y, z = self.read_xyz(basename + ".xyz")
        natoms = len(atoms)

        margin = float(margin)
        extentx = (x.max() - x.min())*margin/2.
        extenty = (y.max() - y.min())*margin/2.
        extentz = (z.max() - z.min())*margin/2.
        margin = int(margin)
        xmin = x.min() * self.ang_to_au - extentx
        xmax = x.max() * self.ang_to_au + extentx
        ymin = y.min() * self.ang_to_au - extenty
        ymax = y.max() * self.ang_to_au + extenty
        zmin = z.min() * self.ang_to_au - extentz
        zmax = z.max() * self.ang_to_au + extentz

        mep_inp = open(basename + "-mep.inp", "w")
        mep_inp.write("{0:d}\n".format(npoints**3))
        for ix in np.linspace(xmin, xmax, npoints, True):
            for iy in np.linspace(ymin, ymax, npoints, True):
                for iz in np.linspace(zmin, zmax, npoints, True):
                    mep_inp.write(
                            "{0:12.6f} {1:12.6f} {2:12.6f}\n".format(
                                ix, iy, iz))
        mep_inp.close()

        orcavpottmp = subprocess.check_output(["orca_vpot", 
                                               basename + ".gbw", 
                                               basename + ".scfp",
                                               basename + "-mep.inp",
                                               basename + "-mep.out"])
        log.debug('Orca_vpot output:\n{}'.format(orcavpottmp))

        vpot = self.read_vpot(basename + "-mep.out")

        #cube = open(basename + "-" + str(npoints/margin) + "-" + str(margin) + "-mep.cube", "w")
        cube = open('{}-{}-{}-mep.cube'.format(basename,npoints/margin,margin), "w")
        cube.write("Generated with ORCA\n")
        cube.write("Electrostatic potential for " + basename + "\n")
        cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
            len(atoms), xmin, ymin, zmin))
        cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
            npoints, (xmax - xmin) / float(npoints - 1), 0.0, 0.0))
        cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
            npoints, 0.0, (ymax - ymin) / float(npoints - 1), 0.0))
        cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}\n".format(
            npoints, 0.0, 0.0, (zmax - zmin) / float(npoints - 1)))
        for i, atom in enumerate(atoms):
            index = self.elements.index(atom)
            cube.write("{0:5d}{1:12.6f}{2:12.6f}{3:12.6f}{4:12.6f}\n".format(
                index, 0.0, x[i] * self.ang_to_au, y[i] * self.ang_to_au, z[i] * self.ang_to_au))

        m = 0 
        n = 0
        vpot = np.reshape(vpot, (npoints, npoints, npoints))
        for ix in range(npoints):
            for iy in range(npoints):
                for iz in range(npoints):
                    cube.write("{0:14.5e}".format(vpot[ix][iy][iz]))
                    m += 1
                    n += 1
                    if (n > 5):
                        cube.write("\n")
                        n = 0
                if n != 0:
                    cube.write("\n")
                    n = 0
        cube.close()
        os.remove(basename + "-mep.inp")
        os.remove(basename + "-mep.out")
        #os.remove(basename + ".K.tmp")
        map(os.remove, glob.glob('*.tmp'))

    def _denscubeinp(self, basename, npoints=80, margin=2):
        '''
        Create a RHO .cube file with ORCA.

        @ basename - file name without the extension;
                     this should be the same for the .gbw and .scfp.
        @ npoints  - number of grid points per original cell size
                     (80 should be fine)

        Dependencies: .xyz .gbw .scfp
        orca has to be properly added to the PATH.
        '''

        npoints = npoints*margin

        atoms, x, y, z = self.read_xyz(basename + ".xyz")
        natoms = len(atoms)

        #extent = margin
        margin = float(margin)
        extentx = (x.max() - x.min())*margin/2.
        extenty = (y.max() - y.min())*margin/2.
        extentz = (z.max() - z.min())*margin/2.
        margin = int(margin)
        xmin = x.min() * self.ang_to_au - extentx
        xmax = x.max() * self.ang_to_au + extentx
        ymin = y.min() * self.ang_to_au - extenty
        ymax = y.max() * self.ang_to_au + extenty
        zmin = z.min() * self.ang_to_au - extentz
        zmax = z.max() * self.ang_to_au + extentz

        densdict = {\
                'npoints': npoints,
                'xmin':    xmin,
                'ymin':    ymin,
                'zmin':    zmin,
                'xmax':    xmax,
                'ymax':    ymax,
                'zmax':    zmax,
                'basename': basename}

        densinp = '''\
! B3LYP cc-pvtz KeepDens #Opt 

%plots
 dim1 {npoints}    # resolution in x-direction
 dim2 {npoints}    # resolution in y-direction
 dim3 {npoints}    # resolution in z-direction
 min1 {xmin}  # x-min value in bohr
 max1 {xmax}  # x-min value in bohr
 min2 {ymin}  # y-min value in bohr
 max2 {ymax}  # y-max value in bohr
 min3 {zmin}  # z-min value in bohr
 max3 {zmax}  # z-max value in bohr
 #
 Format Origin  # Gaussian-cube format 
                       # (an ASCII file)
 ElDens("{basename}-rho.plt");      # Electron density
end

*xyzfile 0 1 {basename}.xyz
'''
        
        rho_inp = open(basename + "-rho.inp", "w")
        rho_inp.write(densinp.format(**densdict))

        #proc = subprocess.Popen(["orca" , "nitromethane-rho.inp"],
        #        stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        #os.system('orca ' + basename + '-rho.inp > '+ basename + '-rho.out')
        #os.system('orca ' + basename + '-rho.inp')
        #subprocess.check_output(['orca' , basename + '-rho.inp'])
        #subprocess.call('orca ' + basename + '-rho.inp > '+ basename + '-rho.out', shell=True)
        #print proc.communicate()
        #rho_out = subprocess.check_output(["orca", basename + "-rho.inp"])
        #with open(basename + "-rho.out", 'wb') as rho_outf:
        #    rho_outf.write(rho_out)

    def cubesize(self, basename, npoints=80, margin=2):

        npoints = npoints*margin

        atoms, x, y, z = self.read_xyz(basename + ".xyz")
        
        natoms = len(atoms)

        #extent = margin
        margin = float(margin)
        extentx = (x.max() - x.min())*margin/2.
        extenty = (y.max() - y.min())*margin/2.
        extentz = (z.max() - z.min())*margin/2.
        margin = int(margin)
        xmin = x.min() * self.ang_to_au - extentx
        xmax = x.max() * self.ang_to_au + extentx
        ymin = y.min() * self.ang_to_au - extenty
        ymax = y.max() * self.ang_to_au + extenty
        zmin = z.min() * self.ang_to_au - extentz
        zmax = z.max() * self.ang_to_au + extentz

        densdict = {\
                'npoints': npoints,
                'xmin':    xmin,
                'ymin':    ymin,
                'zmin':    zmin,
                'xmax':    xmax,
                'ymax':    ymax,
                'zmax':    zmax,
                'margin': margin,
                'np_margin': str(npoints/margin),
                'basename': basename}

        return densdict

    def denscubeinp(self, basename, npoints=80, margin=2):
        '''
        Create a orca_plot input to generate a RHO .cube file with ORCA.

        @ basename - file name without the extension;
                     this should be the same for the .gbw and .scfp.
        @ npoints  - number of grid points per original cell size
                     (80 should be fine)

        '''
        densdict = self.cubesize(basename, npoints, margin)


        densinp = '''\
  2     # PlotType - Type of plot to be done
  7     # Format   - File format for output file
  0   0 # MO and operator (if not density)
  0     # State density to be plotted
{basename}.scfp                              # Input file
{basename}-{np_margin}-{margin}-rho.cube     # Output file
200     # ncont    - number of contours
  1     # icont    - contour option
  0     # Skeleton - flag for Skeleton plotting
  0     # Atoms    - flag for atom plotting
  0     # UseCol   - flag for use of color
  {npoints} {npoints} {npoints} # Grid dimensions
  {xmin}    {xmax}              # X dimension
  {ymin}    {ymax}              # Y dimension
  {zmin}    {zmax}              # Z dimension
 -1  -1  -1                            # defining atoms
    0.000000     0.000000     0.000000 # defining vector-1
    1.000000     0.000000     0.000000 # defining vector-2
    0.000000     1.000000     0.000000 # defining vector-3
'''
        rho_inp = open(basename + "-rho.plot.tmp", "w")
        rho_inp.write(densinp.format(**densdict))

        #proc = subprocess.Popen(["orca" , "nitromethane-rho.inp"],
        #        stderr=subprocess.PIPE, stdout=subprocess.PIPE)
        #os.system('orca ' + basename + '-rho.inp > '+ basename + '-rho.out')
        #os.system('orca ' + basename + '-rho.inp')
        #subprocess.check_output(['orca_plot' , basename + '.gbw', basename + '-rho.plot.tmp'])
        #command = 'orca_plot ' + basename + '.gbw ' + basename + '-rho.plot.tmp'
        #print command
        #os.system(command)
        #subprocess.call('orca ' + basename + '-rho.inp > '+ basename + '-rho.out', shell=True)
        #print proc.communicate()
        #rho_out = subprocess.check_output(["orca", basename + "-rho.inp"])
        #with open(basename + "-rho.out", 'wb') as rho_outf:
        #    rho_outf.write(rho_out)

        
    def _denscube(self, basename, npoints=80, margin=2):
        '''
        Create a RHO .cube file with ORCA.

        @ basename - file name without the extension;
                     this should be the same for the .gbw and .scfp.
        @ npoints  - number of grid points per original cell size
                     (80 should be fine)

        Dependencies: .xyz .gbw .scfp
        orca_plot has to be properly added to the PATH.
        '''
        
        basename = 'nitromethane'
        #self.denscubeinp(basename, npoints, margin)
        #subprocess.check_call(['orca_plot' , basename + '.gbw', basename + '-rho.plot.tmp'])
        #map(os.remove, glob.glob('*.tmp'))

        self._denscubeinp(basename, npoints, margin)
        rho_out = subprocess.check_output(["orca", basename + "-rho.inp"])
        with open(basename + "-rho.out", 'wb') as rho_outf:
            rho_outf.write(rho_out)

    def rhocube(self, basename, npoints=80, margin=2):

        rhodict = self.cubesize(basename, npoints, margin)

        extentx = rhodict['xmax'] - rhodict['xmin']
        extenty = rhodict['ymax'] - rhodict['ymin']
        extentz = rhodict['zmax'] - rhodict['zmin']

        rhodict['sextentx'] = extentx/2.
        rhodict['sextenty'] = extenty/2.
        rhodict['sextentz'] = extentz/2.

        self.ang_to_au = 1.0 / 0.5291772083
        rhodict['xm'] = (rhodict['xmin'] + extentx/2.)/self.ang_to_au
        rhodict['ym'] = (rhodict['ymin'] + extenty/2.)/self.ang_to_au
        rhodict['zm'] = (rhodict['zmin'] + extentz/2.)/self.ang_to_au

        with open(basename + '.inmultiwfn.tmp', 'wb') as multiwfninf:
            multiwfninf.write(multiwfnin.format(**rhodict))

        #os.system('Multiwfn ' + basename + '.wfn < ' + basename + '.inmultiwfn')
        #subprocess.call(['Multiwfn ' , basename + '.wfn < ' + basename + '.inmultiwfn'])
        proc = subprocess.Popen(['Multiwfn', basename + '.wfn'],
                stdin=open(basename + '.inmultiwfn.tmp', 'rb'),
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        proc.communicate()
        #os.rename('density.cub', basename + '-' + str(rhodict['np_margin']) + '-' + str(rhodict['margin']) + '-rho.cube')
        os.rename('density.cub', 
                  '{}-{}-{}-rho.cube'.format(
                  basename,rhodict['np_margin'],rhodict['margin']))
        map(os.remove, glob.glob('*.tmp'))


