#! /usr/bin/octave -q

addpath("../../../src");

mol = mol_readxyz("fragment_1a.xyz");
mol = mol_merge(mol,mol_readxyz("fragment_1b.xyz"));
rep = mol_ball(mol);
rep = mol_stick(mol,rep);
r = [
     0.7970973253 -0.5935178995 0.1112325788 
     0.3703757226 0.3350512385 -0.8663507104 
     0.4769255519 0.7317637801 0.4868927598 
     0.0000000000 0.0000000000 -30.0000000000 
     ];
rep = rep_setdefaultscene(rep,r);
rep_write_obj(rep,"thymine_1.obj");
rep_write_pov(rep,"thymine_1.pov");
system("povray -D -UV +Ithymine_1.pov +Othymine_1.png +W1000 +H1000 +A");

mol = mol_readxyz("fragment_2a.xyz");
mol = mol_merge(mol,mol_readxyz("fragment_2b.xyz"));
mol = mol_merge(mol,mol_readxyz("fragment_2c.xyz"));
rep = mol_ball(mol);
rep = mol_stick(mol,rep);
r = [
     0.8683175445 0.0915508494 -0.4874865115 
     -0.0053238273 0.9844820499 0.1754046381 
     0.4959801733 -0.1497116387 0.8553304076 
     0.0000000000 0.0000000000 -30.0000000000 
     ];
rep = rep_setdefaultscene(rep,r);
rep_write_obj(rep,"thymine_2.obj");
rep_write_pov(rep,"thymine_2.pov");
system("povray -D -UV +Ithymine_2.pov +Othymine_2.png +W1000 +H1000 +A");

