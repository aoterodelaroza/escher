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

class Representation(repi,camangle=np.matrix([80,75,45]),zoom=3,LOG=0):

    '''
    function rep = representation()

    representation - create an empty rep structure.

    Output:
    {rep}: the empty representation.
    '''

    def __init__(self):

        self.name = ""

        self.nball = 0
        self.ball = []

        self.nstick = 0
        self.stick = []

        self.ntriangle = 0
        self.nvertex = 0
        self.triangle = []
        self.vertex = []

        self.cam = type('Cam', (), {})
        self.nlight = 0
        self.light = []
        self.bgcolor = np.zeros(1,3)
    
