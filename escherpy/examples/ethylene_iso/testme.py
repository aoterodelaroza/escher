#!/usr/bin/env python
# coding: utf-8


import escherpy as esc


mol = esc.Molecule()
mol.structfile = esc.escher_data + 'mol/ethylene_iso/c2h4.cube'
isovalue = 0.1
mol.isosurface(esc.escher_data + 'mol/ethylene_iso/c2h4.cube', isovalue)

mol.show()

