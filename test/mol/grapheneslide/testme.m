#! /usr/bin/octave -q

addpath("../../../src");

mol = mol_readxyz("large_00_00.xyz",0);
rep = mol_ball(mol,:,"H",:,:,[200 200 200]);
rep = mol_ball(mol,rep,"C",:,:,[10 10 10]);
rep = mol_stick(mol,rep,"C","H");
rep = mol_stick(mol,rep,"C","C");
rep_write_obj(rep,"large.obj");
r = [
  -0.7160016298 -0.1336150169 0.6851763725 
  0.6641882658 -0.4325072765 0.6097317934 
  0.2148802578 0.8916647434 0.3984378576 
  0.0000000000 0.0000000000 -30.0000000000
];
rep = rep_setdefaultscene(rep,r,0);
rep_write_pov(rep,"large.pov");
system("povray -D -UV +Ilarge.pov +Olarge.png +W1000 +H1000 +A");
