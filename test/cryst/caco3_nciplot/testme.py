#!/usr/bin/env python
# coding: utf-8


import escherpy as esc


mol = esc.Molecule()

mol.structfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
mol.densfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-dens.cube'
mol.gradfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'

mol.readstruct()
mol.stickball()
mol.nciplot()
mol.start()

