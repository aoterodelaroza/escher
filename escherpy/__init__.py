
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

if os.environ['ESCHER_DEBUG']:
    if os.environ['ESCHER_DEBUG'] == 'ON':
        log.setLevel(logging.DEBUG)
    elif os.environ['ESCHER_DEBUG'] == 'OFF':
        log.setLevel(logging.INFO)
else:
    log.setLevel(logging.INFO)


escher_data = os.environ['ESCHER_DATA']

#parser = SafeConfigParser()
#parser.read('escher.ini')

#print parser.get('molecule', 'structfile')
#print parser.get('molecule', 'densfile')
#print parser.get('molecule', 'gradfile')
