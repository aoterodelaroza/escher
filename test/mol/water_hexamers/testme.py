#!/usr/bin/env python
# coding: utf-8


import escherpy as esc


mol = esc.Molecule()
mol.structfile = esc.escher_data + 'mol/water_hexamers/bag.xyz'

mol.readstruct()
#mol.rotx(60.)
#mol.roty(60.)
#mol.rotz(60.)
mol.rot(60., 30., 20.)
mol.stickball()
mol.show()

