#! /etc/alternatives/octave -q

addpath("../../../src");

cr = cr_read_espresso("li.scf.out",0);
mol = cr_crystalbox(cr,[-0.05 -0.25 -0.05],[1.05 1.25 1.05]);
rep = representation();
rep = mol_ball(mol,rep);
rep = mol_stick(mol,rep);
rep = cr_unitcell(cr,rep,:,:,0.03,[0 0 255]);
rep = rep_setdefaultscene(rep);
rep_write_obj(rep,"li.obj");
rep_write_pov(rep,"li.pov");
system("povray -D -UV +Ili.pov +Oli.png +W1000 +H1000 +A");

cr = cr_read_espresso("al.scf.out",0);
mol = cr_crystalbox(cr,[-0.05 -0.25 -0.05],[1.05 1.25 1.05]);
rep = representation();
rep = mol_ball(mol,rep);
rep = mol_stick(mol,rep);
rep = cr_unitcell(cr,rep,0.03,[0 0 255]);
rep = rep_setdefaultscene(rep);
rep_write_obj(rep,"al.obj");
rep_write_pov(rep,"al.pov");
system("povray -D -UV +Ial.pov +Oal.png +W1000 +H1000 +A");
