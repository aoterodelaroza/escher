#!/usr/bin/env python

import os
from math import sqrt
import numpy as np

pwscf_input = \
'''
&control
 calculation='{calculation}',
 restart_mode='from_scratch', 
 title='100', 
 prefix='100',
 pseudo_dir='/home/daniel/pseudo/', 
 outdir=${SCRATCH},
 tprnfor=.true.,
 etot_conv_thr=1e-6, 
 forc_conv_thr=1e-5,
! tstress=.true., 
/
&system
 ibrav=0, 
! celldm(1)={celldm}, 
 nat={nat},
 ntyp=1, 
 ecutwfc={ecutwfc}, 
 ecutrho={ecutrho},
 occupations='smearing',
 smearing='gauss',
 degauss={degauss}
/
&electrons
 conv_thr = 1d-8,
 electron_maxstep = 130,
 mixing_mode  = 'local-TF',
 mixing_beta = 0.3,
 mixing_ndim = 10,
/
&ions
 ion_dynamics='bfgs',
/
&cell
/
CELL_PARAMETERS cubic
 {celldm} 0.000000  0.000000
 0.000000  {celldm2} 0.000000
 0.000000  0.000000  {celldm3}
ATOMIC_SPECIES
 Pt 195.084 {pseudo}
ATOMIC_POSITIONS bohr {positions}
K_POINTS automatic
{n} {n} 3 {s} {s} {s}
'''

qedict = {
        'calculation': 'relax',
	'SCRATCH': '{SCRATCH}',
        'ecutwfc': 120. , # [40, 180]
        'ecutrho': 401.0 ,# [401., 500]
        'degauss': 0.01  ,# [0.01,0.09]
        'n': 6           ,# [4,9]
        's': 1           ,# {0,1}
        'celldm': 5.3273,
        'celldm2': 5.3273,
        'celldm3': 7.5340 ,   # [7.0,9.]
        'pseudo': "Pt.pbe-n-kjpaw_psl.0.1.UPF", # {Pt.pz-n-kjpaw_psl.0.1.UPF, Pt.pbe-n-kjpaw_psl.0.1.UPF}
        'nat': '',
        'positions': ''
}

if os.path.isfile('runme.sh'):
    os.remove('runme.sh')

lldict = {
        'job_type': 'mpich',
        'class': 'plarge',
        'tasks_per_node' : 8,
        'comment': 'details of calculation'
        }
cmd = \
'''#! /bin/bash

#@ job_type       = {job_type}
#@ output         = ll.$(jobid).out
#@ error          = ll.$(jobid).err
#@ class          = {class}
#@ environment    = COPY_ALL
#@ node           = 1
#@ tasks_per_node = {tasks_per_node:d}
#@ comment        = '{comment}'
#@ queue

#echo "HOST: malta"
#echo "SCRATCH dir: $SCRATCH"
#echo "SCRATCH content:"
#ls $SCRATCH

export ESPRESSO_TMPDIR=$SCRATCH

'''

lldict['comment'] = 'vacuum size with 2 cell layers'
with open('runme.sh', 'a') as jobscript:
    jobscript.write(cmd.format(**lldict))

ap = 7.5340
vac = 16.0000
# fcc111
a = ap/sqrt(2.)
b = a
c = ap*sqrt(3.)
pts = [[0., 0., 0.],
       [a/3., 2.*b/3., c/3.],
       [2.*a/3., 1.*b/3., 2.*c/3.]]
qedict['celldm'] = a
qedict['celldm2'] = b
qedict['celldm3'] = c
       
pts = np.array(pts)

for nlayer in range(3):
    pt = pts[nlayer]
    qedict['positions'] = '\n'.join([qedict['positions'],
        'Pt {:.^6.5f} {:.^6.5f} {:.^6.5f}'.format(*pt)])

for nlayer in range(4,16):

    #if nlayer % 2 == 0:
    #    pt = pts[1]
    #elif nlayer % 2 == 1:
    #    pt = pts[0]
    pt = pts[(nlayer-1) % 3]

    pt[2] = pt[2] + c

    qedict['positions'] = '\n'.join([qedict['positions'],
        'Pt {:.^6.5f} {:.^6.5f} {:.^6.5f}'.format(*pt)])

    print qedict['positions']

    qedict['celldm3'] = float(vac) + float(nlayer-1)*0.3333333*c
    print qedict['celldm3'] 

    qedict['nat'] = nlayer

    with open('{}.in'.format(nlayer), 'w') as pwscfin:
        pwscfin.write(pwscf_input.format(**qedict))

for nlayer in range(4,16):

    with open('runme.sh', 'a') as jobscript:
        jobscript.write('mpirun -np {0} pw.x < {1}.in > {1}.out \n'.format(
                                        lldict['tasks_per_node'], nlayer))




os.system('chmod u+x runme.sh')

if os.path.isfile('llsend'):
    os.remove('llsend')

llsend = \
'''#! /bin/bash

echo "vacuum analysis" >> sent
llsubmit runme.sh >> sent

'''

with open('llsend', 'a') as jobscript:
    jobscript.write(llsend)

os.system('chmod u+x llsend')
