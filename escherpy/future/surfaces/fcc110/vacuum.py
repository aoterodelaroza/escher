#!/usr/bin/env python

import os

pwscf_input = \
'''
&control
 restart_mode='from_scratch', 
 title='100', 
 prefix='100',
 pseudo_dir='/home/daniel/pseudo/', 
 outdir=${SCRATCH},
! tstress=.true., 
 tprnfor=.true.,
 calculation='scf',
 etot_conv_thr=1e-6, 
 forc_conv_thr=1e-5,
/
&system
 ibrav=0, 
! celldm(1)={celldm}, 
 nat=5,
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
/
&cell
/
CELL_PARAMETERS cubic
 {celldm} 0.000000  0.000000
 0.000000  {celldm2} 0.000000
 0.000000  0.000000  {celldm3}
ATOMIC_SPECIES
 Pt 195.084 {pseudo}
ATOMIC_POSITIONS bohr
Pt  0.0000  0.0000  0.0000
Pt  2.6637  2.6637  3.7670
Pt  0.0000  5.3273  0.0000
Pt  2.6637  7.9911  3.7670
Pt  0.0000  10.655  0.0000
K_POINTS automatic
{n} {n} 3 {s} {s} {s}
'''

qedict = {
	'SCRATCH': '{SCRATCH}',
        'ecutwfc': 120. , # [40, 180]
        'ecutrho': 401.0 ,# [401., 500]
        'degauss': 0.01  ,# [0.01,0.09]
        'n': 6           ,# [4,9]
        's': 1           ,# {0,1}
        'celldm': 5.3273,
        'celldm2': 5.3273,
        'celldm3': 7.5340 ,   # [7.0,9.]
        'pseudo': "Pt.pbe-n-kjpaw_psl.0.1.UPF" # {Pt.pz-n-kjpaw_psl.0.1.UPF, Pt.pbe-n-kjpaw_psl.0.1.UPF}
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

'''

lldict['comment'] = 'vacuum size with 2 cell layers'
with open('runme.sh', 'a') as jobscript:
    jobscript.write(cmd.format(**lldict))


cbulk = qedict['celldm2']
for vac in range(10,19,2):
    qedict['celldm2'] = float(vac) + 2*cbulk

    with open('{}.in'.format(vac), 'w') as pwscfin:
        pwscfin.write(pwscf_input.format(**qedict))

    with open('runme.sh', 'a') as jobscript:
        jobscript.write('mpirun -np {0} pw.x < {1}.in > {1}.out \n'.format(lldict['tasks_per_node'], vac))

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
