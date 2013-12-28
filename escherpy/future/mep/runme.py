#!/usr/bin/env python

from mep import MEP

mep = MEP()

for margin in range(1,4):
    for np in range(60, 100, 20):
        mep.orca2potcube('nitromethane', np, margin)
