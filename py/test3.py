#!/usr/bin/env python
# coding: utf-8

import os
from molecule import Molecule

escher_data = os.environ['ESCHER_DATA']

mol3 = Molecule()

mol3.structfile = escher_data + 'mol/water_hexamers/bag.xyz'
#mol3.structfile =  '/home/daniel/calc/espresso/pt/fcc_100.xyz'
mol3.readstruct()
mol3.stickball()

mol3.start()

