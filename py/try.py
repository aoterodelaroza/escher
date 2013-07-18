# coding: utf-8
import os
from molecule import Molecule

escher_data = os.environ['ESCHER_DATA']

mol = Molecule()
#mol.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
#a.readcube('/home/daniel/pkg/vis/VTKData/Data/m4_TotalDensity.cube')
mol.readcube(escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube')
mol.stickball()
#mol.isosurface('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-grad.cube', 0.5)
mol.nciplot(escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube',
            escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube')
mol.start()

#mol.readcube('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-grad.cube')
#mol.nciplot('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-dens.cube',
#            '/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-grad.cube')
