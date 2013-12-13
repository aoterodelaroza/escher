#!/usr/bin/env python
# coding: utf-8


import escherpy as esc


mol = esc.Molecule()

mol.structfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
mol.readstruct()
mol.stickball()

densfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-dens.cube'
gradfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
mol.nciplot(densfile, gradfile)

mol.show()

