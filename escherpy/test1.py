#!/usr/bin/env python
# coding: utf-8

import os
from molecule import Molecule

escher_data = os.environ['ESCHER_DATA']

mol = Molecule()
mol.structfile = escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
mol.densfile = escher_data + 'cryst/aragonite/aragonite.9.4-dens.cube'
mol.gradfile = escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'

mol.readstruct()
mol.stickball()
mol.nciplot()


mol.start()


#mol.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
#a.readcube('/home/daniel/pkg/vis/VTKData/Data/m4_TotalDensity.cube')
#mol.readcube('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-grad.cube')
#mol.nciplot('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-dens.cube',
#            '/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-grad.cube')
