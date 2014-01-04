#!/usr/bin/env python

import os
import subprocess
import numpy as np
from templates.critic2 import rhoinp, mepinp


class Critic2(object):


    def __init__(self):
        pass



    def _rho(self, critic2dict):


        basename = critic2dict['basename']
        np_margin = critic2dict['np_margin']
        margin = critic2dict['margin']

        rhoinpf = open('{}-{}-{}-rho.incritic'.format(
                       basename,np_margin,margin), "w")
        rhoinpf.write(rhoinp.format(**critic2dict))


    def rho(self, basename, np_margin, margin):

        critic2dict = {\
                'npoints'   : np_margin*margin,
                'margin'    : margin,
                'np_margin' : np_margin,
                'basename'  : basename}

        np_margin = critic2dict['np_margin']
        margin = critic2dict['margin']

        self._rho(critic2dict)
        critic2out = subprocess.check_output(
                    ['critic2' , 
                    '{}-{}-{}-rho.incritic'.format(
                    basename,np_margin,margin)])
        with open('{}-{}-{}-rho.outcritic'.format(
                  basename,np_margin,margin), 
                  'wb') as critic2outf:
            critic2outf.write(critic2out)


    def _mep(self, critic2dict):


        basename = critic2dict['basename']
        np_margin = critic2dict['np_margin']
        margin = critic2dict['margin']

        mepinpf = open('{}-{}-{}-mep.incritic'.format(
                       basename,np_margin,margin), "w")
        mepinpf.write(mepinp.format(**critic2dict))


    def mep(self, basename, np_margin, margin):

        critic2dict = {\
                'npoints'   : np_margin*margin,
                'margin'    : margin,
                'np_margin' : np_margin,
                'basename'  : basename}

        np_margin = critic2dict['np_margin']
        margin = critic2dict['margin']

        self._mep(critic2dict)
        critic2out = subprocess.check_output(
                    ['critic2' , 
                    '{}-{}-{}-mep.incritic'.format(
                    basename,np_margin,margin)])
        with open('{}-{}-{}-mep.outcritic'.format(
                  basename,np_margin,margin), 
                  'wb') as critic2outf:
            critic2outf.write(critic2out)

    def openf(self, filename):
        """
        Opens a file and returns lines without whitespaces.
        """
        with open(filename, 'rb') as file:
            text = file.read()

        lines = [line for line in text.splitlines()]
        lines = map(lambda s: s.strip(), lines)
        return lines

    def readoutcritic(self, filename):
        """
        Extracts basin charges, ...
        from a *.outcritic file
        """
        lines = self.openf(filename)

        charges = []
        volumes = []

        for n,i in enumerate(lines):

            if i.startswith('Number of atoms in the unit cell'):
                self.nat = int(i.split()[-1])

            if i.startswith('* List of basins and local properties'):
                for j in range(2,self.nat+1):
                    charges.append(float(lines[j+n].split()[-1]))
                    volumes.append(float(lines[j+n].split()[-2]))
                self.charges = charges
                self.volume = volumes

            if i.startswith('Total integration'):
                self.totcharge = float(i.split()[-1])


