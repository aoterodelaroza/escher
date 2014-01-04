#!/usr/bin/env python

import shutil
import matplotlib.pyplot as plt
from escherpy import escher_data
from escherpy.orca import Orca
from escherpy.critic2 import Critic2    

orca = Orca()
critic2 = Critic2()

basename = 'nitromethane'
shutil.copy2('{0}mol/{1}/{1}.xyz'.format(escher_data,basename), 
             '{}.xyz'.format(basename))

orca.run(basename)
orca.genwfn(basename)

nestedcharges = []
totcharges = []
for margin in range(1,3):
    for np in range(10, 20, 20):
        orca.potcube(basename, np, margin)
        orca.rhocube(basename, np, margin)
        critic2.rho(basename, np, margin)
        critic2.mep(basename, np, margin)
        critic2.readoutcritic('{}-{}-{}-mep.outcritic'.format(basename,np,margin))
        nestedcharges.append(critic2.charges)
        totcharges.append(critic2.totcharge)
        pass


def basinpropplot():

    x = range(1, 3)
    charges = zip(*nestedcharges)
    for i in charges:
        plt.plot(x,i)

    plt.show()


    plt.plot(x,totcharges)
    plt.show()


basinpropplot()




#subprocess.check_call(['orca_plot' , basename + '.gbw', basename + '-rho.plot.tmp'])


#rho_out = subprocess.check_output(['orca', 'nitromethane-rho.inp'])
#with open("nitromethane-rho.out", 'wb') as rho_outf:
#    rho_outf.write(rho_out)
