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

class Grid():

    '''
    function g = grid()

    grid - create an empty grid structure.

    Output:
    {g}: the empty grid structure with all the fields defined.
    '''

    def __init__(self):

        self.name = ""
        self.x0 = np.matrix([0,0,0])
        self.dx = np.zeros(3,3)
        self.a = np.zeros(3,3)
        self.n = np.matrix([0,0,0])
        self.f = []
        self.omega = 0.0

