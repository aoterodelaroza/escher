#! /etc/alternatives/octave -q

addpath("../../../src");
cr = cr_read_vasp("POSCAR","POTCAR");
cr_write_tessel(cr,"doped_si.tess");
cr_write_critic2(cr,"doped_si.incritic");

## crystalbox but remove the isolated atoms
mol = cr_crystalbox(cr,:,:,2); 

## Re-center molecule
cm = mol_cmass(mol,:,0);
for i = 1:length(mol.atname)
  mol.atxyz(:,i) -= cm;
endfor

## make the plot
rep = mol_ball(mol);
rep = mol_stick(mol,rep);
r = [
     0.9784541726 -0.0114633143 0.2061457932 -100.0078125000
     0.2033965737 0.2249962687 -0.9528941512 -100.0000000000
     -0.0354587957 0.9742925763 0.2224801183 -100.0000000000
     0.0000000000 0.0000000000 -30.0000000000 0.0000000000
     ];
rep = rep_setdefaultscene(rep,r);
rep_write_obj(rep,"siph_a.obj");
rep_write_pov(rep,"siph_a.pov");
system("povray -D -UV +Isiph_a.pov +Osiph_a.png +W1000 +H1000 +A");


