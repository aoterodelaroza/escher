#!/usr/bin/env python
# coding: utf-8


import escherpy as esc


mol = esc.Molecule()

mol.structfile = esc.escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'

mol.basinfile = esc.escher_data + 'cryst/aragonite/aragonite-1.basin'
mol.readsurf('Ca')
mol.basinfile = esc.escher_data + 'cryst/aragonite/aragonite-5.basin'
mol.readsurf('C')
mol.basinfile = esc.escher_data + 'cryst/aragonite/aragonite-9.basin'
mol.readsurf('O')
mol.basinfile = esc.escher_data + 'cryst/aragonite/aragonite-13.basin'
mol.readsurf('O')

mol.readstruct()
mol.stickball()
#mol.nciplot()
mol.start()

