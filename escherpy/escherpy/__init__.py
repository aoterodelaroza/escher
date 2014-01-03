
__doc__ = '''\
Escherpy is a python package for tranforming quantum chemistry common files,
manipulate molecular/crystal data and visualize it.
'''
__author__ = 'Daniel Menendez'

# from escherpy import *
# * means __all__ modules
__all__ = ['molecule', 'crystal']

import os
#from ConfigParser import SafeConfigParser
from molecule import Molecule
from crystal import Crystal
import logging
import logcolors

log = logcolors.escherlog()

debug = os.environ.get('ESCHER_DEBUG')
if debug: 
    if debug == 'ON':
        log.setLevel(logging.DEBUG)
    elif debug == 'OFF':
        log.setLevel(logging.INFO)
else:
    log.setLevel(logging.INFO)
del debug

#try: 
#    debug = os.environ['ESCHER_DEBUG']
#    if debug == 'ON':
#        log.setLevel(logging.DEBUG)
#    elif debug == 'OFF':
#        log.setLevel(logging.INFO)
#except KeyError:
#    log.setLevel(logging.INFO)
#    raise KeyError

escher_data = os.environ.get('ESCHER_DATA')
if not escher_data:
    print 'You need to define the shell variable ESCHER_DATA pointing to your escher dir'
    exit()

#try:
#    escher_data = os.environ['ESCHER_DATA']
#except KeyError:
#    print 'You need to define the shell variable ESCHER_DATA pointing to your escher dir'
#    raise KeyError

#parser = SafeConfigParser()
#parser.read('../escher.ini')

#print parser.get('molecule', 'structfile')
#print parser.get('molecule', 'densfile')
#print parser.get('molecule', 'gradfile')
