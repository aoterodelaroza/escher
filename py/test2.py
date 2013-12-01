#!/usr/bin/env python
# coding: utf-8

import os
from molecule import Molecule

escher_data = os.environ['ESCHER_DATA']


mol2 = Molecule()
mol2.structfile = escher_data + 'mol/ethylene_iso/c2h4.cube'
mol2.isovalue = 0.1
mol2.isosurface(escher_data + 'mol/ethylene_iso/c2h4.cube')

mol2.start()

