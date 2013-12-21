# Copyright (C) 2011 Victor Lua~na and Alberto Otero-de-la-Roza
#
# This octave routine is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version. See <http://www.gnu.org/licenses/>.
#
# The routine distributed in the hope that it will be useful, but WITHOUT
# ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
# more details.

from logging import getLogger
import numpy as np
from molecule import Molecule
from grid import Grid

log = getLogger('escherlog')

class Crystal(Molecule):

    '''
    function cr = crystal()

    crystal - create an empty crystal structure and initialize the number
    of atoms to zero.

    Output:
    {cr}: the empty crystal structure with all the fields defined.

    '''

    def __init__(self):

        self.name    = ""
        self.nat     = 0
        self.ntyp    = 0
        self.attyp   = []
        self.rvdwtyp = []
        self.c6typ   = []
        self.zvaltyp = []
        self.ztyp    = []
        self.typ     = []
        self.x       = []
        self.cellparams   = np.zeros(6)
        self.omega   = NaN

        self.grid = Grid()

        super(Crystal, self).__init__()
        # python 3.x
        #super().__init__()

    def cellbox(self):

        stickbond = self.rep.stick(self.atxyz[i], self.atxyz[j])
        self.rep.ren.AddActor(stickbond)

    def steitz_op(self, rot, translation):
        pass

