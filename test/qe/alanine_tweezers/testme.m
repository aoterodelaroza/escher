#! /etc/alternatives/octave -q

addpath("../../../src");

## colors and grayscale definition
rgbc = [10, 10, 10];
rgbh = [200,200,200];
rgbo = [255,013,013];
rgbn = [048,080,248];
function x = gray(r)
  x = 0.21 * r(1) + 0.71 * r(2) + 0.07 * r(3);
endfunction

## racemate
cr = cr_read_espresso("ala_dl.scf.out",0);
## green -> L 
mol1 = mol_readxyz("ala_dl-fragment_l.xyz");
rep = representation();
rep = mol_ball(mol1,rep,"C",0,:,[0 gray(rgbc) 0]);
rep = mol_ball(mol1,rep,"H",0,:,[0 gray(rgbh) 0]);
rep = mol_ball(mol1,rep,"O",0,:,[0 gray(rgbo) 0]);
rep = mol_ball(mol1,rep,"N",0,:,[0 gray(rgbn) 0]);
## red -> D
mol2 = mol_readxyz("ala_dl-fragment_d.xyz");
rep = mol_ball(mol2,rep,"C",0,:,[gray(rgbc) 0 0]);
rep = mol_ball(mol2,rep,"H",0,:,[gray(rgbh) 0 0]);
rep = mol_ball(mol2,rep,"O",0,:,[gray(rgbo) 0 0]);
rep = mol_ball(mol2,rep,"N",0,:,[gray(rgbn) 0 0]);
## plot all
mol = mol_merge(mol1,mol2);
rep = mol_stick(mol,rep);
rep = mol_stick(mol,rep,"H","O",[1.6 2.0],:,0.01,[0 200 0]);
rep = cr_unitcell(cr,rep,:,:,0.03,[0 0 255]);
r = [
     -0.0144468546 -0.9801896214 -0.1975281090 0.0
     -0.9997020364 0.0102860928 0.0220856741   0.0
     -0.0196165890 0.1977885664 -0.9800460339  0.0
     0.0000000000 0.0000000000 -33.0000000000  0.0
     ];
rep = rep_setdefaultscene(rep,r);
rep_write_obj(rep,"ala_dl.obj");
rep_write_pov(rep,"ala_dl.pov");
system("povray -D -UV +Iala_dl.pov +Oala_dl.png +W1000 +H1000 +A");

## conglomerate
cr = cr_read_espresso("ala_l.scf.out",0);
mol = mol_readxyz("ala_l.xyz");
rep = representation();
rep = mol_ball(mol,rep,"C",0,:,[0 gray(rgbc) 0]);
rep = mol_ball(mol,rep,"H",0,:,[0 gray(rgbh) 0]);
rep = mol_ball(mol,rep,"O",0,:,[0 gray(rgbo) 0]);
rep = mol_ball(mol,rep,"N",0,:,[0 gray(rgbn) 0]);
rep = mol_stick(mol,rep);
rep = mol_stick(mol,rep,"H","O",[1.6 2.0],:,0.01,[0 200 0]);
rep = cr_unitcell(cr,rep,0.03,[0 0 255]);
r = [
     0.9999222159 -0.0062757158 0.0107810823 
     0.0015675970 0.9206032753 0.3904960155 
     -0.0123757413 -0.3904487491 0.9205414653
     0.0000000000 0.0000000000 -40.0000000000
     ];
rep = rep_setdefaultscene(rep,r);
rep_write_obj(rep,"ala_l.obj");
rep_write_pov(rep,"ala_l.pov");
system("povray -D -UV +Iala_l.pov +Oala_l.png +W1000 +H1000 +A");
