#!/usr/bin/env python

import os
import subprocess
import numpy as np


class MEP(object):

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

    def read_xyz(self, xyz):
        atoms = []
        x = y = z = []
        f = open(xyz, "r")
        f.next()
        f.next()

        for line in f:
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

    def orca2potcube(self, basename, npoints=80, margin=2):
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

        mep_inp = open(basename + "-mep.inp", "w")
        mep_inp.write("{0:d}\n".format(npoints**3))
        for ix in np.linspace(xmin, xmax, npoints, True):
            for iy in np.linspace(ymin, ymax, npoints, True):
                for iz in np.linspace(zmin, zmax, npoints, True):
                    mep_inp.write("{0:12.6f} {1:12.6f} {2:12.6f}\n".format(ix, iy, iz))
        mep_inp.close()

        subprocess.check_call(["orca_vpot", basename + ".gbw", basename + ".scfp",
                basename + "-mep.inp", basename + "-mep.out"])

        vpot = self.read_vpot(basename + "-mep.out")

        cube = open(basename + "-" + str(npoints/margin) + "-" + str(margin) + "-mep.cube", "w")
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

        m = n = 0
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
        os.remove(basename + ".K.tmp")

