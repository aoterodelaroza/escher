#! /usr/bin/octave -q

addpath("../../../src");

cr = cr_read_espresso("BZDMAZ01_withH.scf.out",0);
mol = cr_molmotif(cr);
rep = representation();
rep = mol_ball(mol,rep);
rep = mol_stick(mol,rep);
rep = cr_unitcell(cr,rep,:,:,0.03,[0 0 255]);
r = [
     1.0000000000 0.0000000000 0.0000000000 -2.1762695345
     0.0000000000 1.0000000000 0.0000000000 0.0000000000
     0.0000000000 0.0000000000 1.0000000000 0.0000000000
     0.0000000000 0.0000000000 -60.0000000000 0.0000000000
     ];
rep = rep_setdefaultscene(rep,r,0);
rep_write_obj(rep,"BZDMAZ01.obj");
rep_write_pov(rep,"BZDMAZ01.pov");
system("povray -D -UV +IBZDMAZ01.pov +OBZDMAZ01.png +W1000 +H1000 +A");
