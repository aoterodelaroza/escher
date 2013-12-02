
__all__ = ['molecule']

import os
#from ConfigParser import SafeConfigParser
from molecule import Molecule
import logcolors

log = logcolors.escherlog()
#log.setLevel(log.DEBUG)

escher_data = os.environ['ESCHER_DATA']

#parser = SafeConfigParser()
#parser.read('escher.ini')

#print parser.get('molecule', 'structfile')
#print parser.get('molecule', 'densfile')
#print parser.get('molecule', 'gradfile')
