#! /usr/bin/octave -q

addpath("../../../src");

r = [
     0.8607525229, 0.2052026689, -0.4658299088,
     0.0513425246, 0.8754761219, 0.4805266559,
     0.5064282417, -0.4375315011, 0.7430326939,
     0.0000000000, 0.0000000000, -15
     ];

atoms = {"al","p","s","si"};
for i = 1:length(atoms)
  at = atoms{i};
  mol = mol_readcube(sprintf("%s-grad.cube",at));
  rdg = grid_readcube(sprintf("%s-grad.cube",at));
  rep = representation();
  rep = mol_ball(mol,rep,"O");
  rep = mol_ball(mol,rep,"H");
  rep = mol_stick(mol,rep);
  rep = grid_isosurface(rdg,rep,0.5,:,:,:,:,0.002);
  rep_write_obj(rep,sprintf("%s.obj",at));
  rep = rep_setdefaultscene(rep,r);
  rep_write_pov(rep,sprintf("%s.pov",at));
  system(sprintf("povray -D -UV +I%s.pov +O%s.png +W2000 +H2000 +A",at,at));
endfor
