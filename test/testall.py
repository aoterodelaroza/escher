#!/usr/bin/env python

'''
This is an script to search and run all
escherpy tests in this directory
'''

import os

pwd = os.getcwd()

#testlist = \
#        ['/cryst/caco3_nciplot/',
#        '/cryst/caco3_basin/',
#        '/mol/ethylene_iso/',
#        '/mol/nitromethane/',
#        '/mol/water_hexamers/']

testlist = []
for root, dirs, files in os.walk(os.getcwd()):
    for file in files:
        if file.endswith('testme.py'):
            testlist.append(root)

for testdir in testlist:
    #os.chdir(os.path.join(os.getcwd() + testdir))
    os.chdir(testdir)
    os.system('./testme.py')
    os.chdir(pwd)

