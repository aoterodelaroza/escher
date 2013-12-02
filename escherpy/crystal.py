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

import numpy as np

class Crystal():

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
        self.a       = np.zeros(1,3)
        self.b       = np.zeros(1,3)
        self.omega   = 0.0

