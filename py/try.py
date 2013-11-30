# coding: utf-8
import os
from molecule import Molecule

escher_data = os.environ['ESCHER_DATA']

#mol = Molecule()
#mol.structfile = escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
#mol.densfile = escher_data + 'cryst/aragonite/aragonite.9.4-dens.cube'
#mol.gradfile = escher_data + 'cryst/aragonite/aragonite.9.4-grad.cube'
#
#mol.readstruct()
#mol.stickball()
#mol.nciplot()
#
#
#mol.start()

# ====================================================================================
#mol2 = Molecule()
#mol2.structfile = escher_data + 'mol/ethylene_iso/c2h4.cube'
#mol2.isovalue = 0.1
#mol2.isosurface(escher_data + 'mol/ethylene_iso/c2h4.cube')

#mol2.start()

# ====================================================================================
mol3 = Molecule()

mol3.structfile = escher_data + 'mol/water_hexamers/bag.xyz'
#mol3.structfile =  '/home/daniel/calc/espresso/pt/fcc_100.xyz'
mol3.readstruct()
mol3.stickball()

mol3.start()

#mol.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
#a.readcube('/home/daniel/pkg/vis/VTKData/Data/m4_TotalDensity.cube')
#mol.readcube('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-grad.cube')
#mol.nciplot('/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-dens.cube',
#            '/home/daniel/calc/espresso/CaCO3/e_v2/arag/aragonite.9.4-grad.cube')
