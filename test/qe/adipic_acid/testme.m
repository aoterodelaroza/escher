#! /etc/alternatives/octave -q

addpath("../../../src");

cr = cr_read_espresso("adipic_acid.scf.out",0);
mol = cr_molmotif(cr);
rep = representation();
rep = mol_ball(mol,rep);
rep = mol_stick(mol,rep);
rep = cr_unitcell(cr,rep,:,:,0.03,[0 0 255]);
rep = rep_setdefaultscene(rep);
rep_write_pov(rep,"adipic_acid_cell.pov");
system("povray -D -UV +Iadipic_acid_cell.pov +Oadipic_acid_cell.png +W1000 +H1000 +A");
