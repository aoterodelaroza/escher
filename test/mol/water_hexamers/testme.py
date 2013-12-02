#!/usr/bin/env python
# coding: utf-8


import escherpy as esc


mol = esc.Molecule()
mol.structfile = esc.escher_data + 'mol/water_hexamers/bag.xyz'
mol.readstruct()
mol.stickball()

mol.start()

