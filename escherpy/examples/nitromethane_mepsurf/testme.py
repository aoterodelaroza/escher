#!/usr/bin/env python
# coding: utf-8

# A MEP isosurface


import escherpy as esc



mol = esc.Molecule()

mol.structfile = esc.escher_data + 'mol/nitromethane/nitromethane.xyz'
mol.readstruct()
mol.stickball()

isovalue = -0.05
mol.isosurface(esc.escher_data + 'mol/nitromethane/potential.cube', isovalue)

#mol.readcps()
#mol.cpball()
#mol.nciplot()
mol.show()
