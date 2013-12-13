#!/usr/bin/env python
# coding: utf-8

# A multiwfn critical points analyzer


import escherpy as esc



mol = esc.Molecule()

mol.structfile = esc.escher_data + 'mol/nitromethane/nitromethane.xyz'
isovalue = -0.05
mol.isosurface(esc.escher_data + 'mol/nitromethane/potential.cube', isovalue)

mol.readstruct()
mol.stickball()
#mol.readcps()
#mol.cpball()
#mol.nciplot()
mol.show()
