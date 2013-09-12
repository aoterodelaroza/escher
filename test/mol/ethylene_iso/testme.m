#! /usr/bin/octave -q

r = [
     0.1135832071, 0.1106972098, -0.9873424768,
     -0.0290711969, 0.9937183857, 0.1080677286,
     0.9931031466, 0.0164285377, 0.1160877943,
     0.0000000000, 0.0000000000, -7.0000000000
     ];

mol = mol_readcube("c2h4.cube");
grid = grid_readcube("c2h4.cube");

rep = mol_ball(mol);
rep = mol_stick(mol,rep);
rep = grid_isosurface(grid,rep,0.1,[255 0 0 225],[128 0 0]);

rep = rep_setdefaultscene(rep,r);
rep_write_obj(rep,"c2h4.obj");
rep_write_pov(rep,"c2h4.pov");
system("povray -D -UV +Ic2h4.pov +Oc2h4.png +W2000 +H2000 +A");
