#!/usr/bin/env python
# coding: utf-8

# A multiwfn critical points analyzer


import escherpy as esc


mol = esc.Molecule()

mol.structfile = esc.escher_data + 'mol/nitromethane/nitromethane.xyz'
mol.readstruct()
mol.stickball()

mol.cpsfile = esc.escher_data + 'mol/nitromethane/CPprop.txt'
mol.readcps()
mol.cpball()

mol.show()

