#!/usr/bin/env python
# coding: utf-8

# A multiwfn critical points analyzer


import escherpy as esc



mol = esc.Molecule()

mol.structfile = esc.escher_data + 'mol/nitromethane/nitromethane.xyz'
mol.readstruct()
mol.stickball()

isovalue = 0.02
isofile = esc.escher_data + 'mol/nitromethane/density.cube'
mapfile = esc.escher_data + 'mol/nitromethane/potential.cube'
rangecolor = [-0.01, 0.01]
colortable = [[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]]
mol.surfmap(mapfile, isofile, isovalue, rangecolor, colortable)

mol.show()
