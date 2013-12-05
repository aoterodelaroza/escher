#!/usr/bin/env python
# coding: utf-8

# A multiwfn critical points analyzer


import escherpy as esc


mol = esc.Molecule()

mol.structfile = 'nitromethane.xyz'
mol.cpsfile = 'CPprop.txt'

mol.readstruct()
mol.stickball()
mol.readcps()
mol.cpball()
#mol.nciplot()
mol.start()

