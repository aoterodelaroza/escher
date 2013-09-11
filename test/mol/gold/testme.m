#! /etc/alternatives/octave -q

addpath("../../../src");

mol = mol_readxyz("ortho.xyz",0);
rep = mol_ball(mol);
rep = mol_stick(mol,rep,"P","H",:,0);
rep = mol_stick(mol,rep,"Au","P",:,0);
rep = mol_stick(mol,rep,"Au","Cl",:,0);
rep_write_obj(rep,"ortho.obj");
r = [
     -0.2186212540 -0.1323426366 0.9667922258 
     0.0243301652 -0.9911892414 -0.1301806569 
     0.9755048752 -0.0049382150 0.2199179530 
     0.0000000000 0.0000000000 -20.0000000000
     ];
rep = rep_setdefaultscene(rep,r,0);
rep_write_pov(rep,"ortho.pov");
system("povray -D -UV +Iortho.pov +Oortho.png +W1000 +H1000 +A");

mol = mol_readxyz("para.xyz",0);
rep = mol_ball(mol);
rep = mol_stick(mol,rep,"P","H",:,0);
rep = mol_stick(mol,rep,"Au","P",:,0);
rep = mol_stick(mol,rep,"Au","Cl",:,0);
rep_write_obj(rep,"para.obj");
rep = rep_setdefaultscene(rep,:,0);
rep_write_pov(rep,"para.pov");
system("povray -D -UV +Ipara.pov +Opara.png +W1000 +H1000 +A");

mol = mol_readxyz("eta2.xyz",0);
rep = mol_ball(mol);
rep = mol_stick(mol,rep,"P","H",:,0);
rep = mol_stick(mol,rep,"Au","P",:,0);
rep = mol_stick(mol,rep,"Au","Cl",:,0);
rep_write_obj(rep,"eta2.obj");
r = [
     0.0 1.0 0.0 
     0.0 0.0 1.0 
     1.0 0.0 0.0 
     0.0 0.0 -20.
     ];
rep = rep_setdefaultscene(rep,r,0);
rep_write_pov(rep,"eta2.pov");
system("povray -D -UV +Ieta2.pov +Oeta2.png +W1000 +H1000 +A");

mol = mol_readxyz("trimer.xyz",0);
rep = mol_ball(mol);
rep = mol_stick(mol,rep,"P","H",:,0);
rep = mol_stick(mol,rep,"Au","P",:,0);
rep = mol_stick(mol,rep,"Au","Cl",:,0);
rep_write_obj(rep,"trimer.obj");
r = [
     -0.0603359938 0.8680983186 0.4927101135 
     -0.4277157187 -0.4684816599 0.7730342746
     0.9018964767 -0.1640976071 0.3995670676 
     0.0000000000 0.0000000000 -25.0000000000
     ];
rep = rep_setdefaultscene(rep,r,0);
rep_write_pov(rep,"trimer.pov");
system("povray -D -UV +Itrimer.pov +Otrimer.png +W1000 +H1000 +A");

mol = mol_readxyz("cube.xyz",0);
rep = mol_ball(mol);
rep = mol_stick(mol,rep,"P","H",:,0);
rep = mol_stick(mol,rep,"Au","P",:,0);
rep = mol_stick(mol,rep,"Au","Cl",:,0);
rep_write_obj(rep,"cube.obj");
r = [
     0.9282323718 -0.1410374045 0.3442282677 
     0.0067285821 0.9315567017 0.3635338545 
     -0.3719400167 -0.3351277113 0.8656500578
     0.0000000000 0.0000000000 -30.0000000000
     ];
rep = rep_setdefaultscene(rep,r,0);
rep_write_pov(rep,"cube.pov");
system("povray -D -UV +Icube.pov +Ocube.png +W1000 +H1000 +A");

cr = cr_read_espresso("aucl_tet.scf.out",0);
mol = cr_crystalbox(cr);
rep = mol_ball(mol);
rep = cr_unitcell(cr,rep,:,:,:,[0 0 255]);
rep = mol_stick(mol,rep,"Au","Cl",:,0);
rep_write_obj(rep,"crystal.obj");
rep = rep_setdefaultscene(rep,:,0);
rep_write_pov(rep,"crystal.pov");
system("povray -D -UV +Icrystal.pov +Ocrystal.png +W1000 +H1000 +A");

