# coding: utf-8
from molecule import Molecule


mol = Molecule()
mol.readcube('../../escher_data/mol/ethylene_iso/c2h4.cube')
#a.readcube('/home/daniel/pkg/vis/VTKData/Data/m4_TotalDensity.cube')
#a.readcube('/home/daniel/calc/espresso/CaCO3/basins/arag/aragonite.8.8.cube')
mol.nciplot('a', 'b')
mol.stickball()
mol.isosurface()
mol.start()
