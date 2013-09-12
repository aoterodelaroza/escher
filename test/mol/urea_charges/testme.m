#! /usr/bin/octave -q

addpath("../../../src");

cr = cr_read_espresso("urea.scf.out");
mol = mol_readxyz("urea.xyz");

function rgb = colorscale(q,m)
  rgb = [255 255 255 0 0];
  if (q > 0)
    rgb(1:2) -= floor((q / m) * 255);
  else
    rgb(2:3) -= floor(-(q / m) * 255);
  endif
endfunction

qs = [0.48309326 1.35048571 -1.17798673 -0.88775211];

m = max(abs(qs));
rep = representation();
rep = mol_ball(mol,rep,"H",:,:,colorscale(qs(1),m));
rep = mol_ball(mol,rep,"C",:,:,colorscale(qs(2),m));
rep = mol_ball(mol,rep,"N",:,:,colorscale(qs(3),m));
rep = mol_ball(mol,rep,"O",:,:,colorscale(qs(4),m));
rep = mol_stick(mol,rep,:,:,:,:,0.07,[0 0 0]);
rep = cr_unitcell(cr,rep,:,:,:,[0 0 0]);
r = [
     0.9886412621 0.0500601456 -0.1417127401
     -0.1502810121 0.3420233727 -0.9275969863
     0.0020335019 0.9383574128 0.3456615210 
     0.0000000000 0.0000000000 -30.0000000000
     ];
rep = rep_setdefaultscene(rep,r);
rep_write_obj(rep,"urea.obj");
rep_write_pov(rep,"urea.pov");
system("povray -D -UV +Iurea.pov +Ourea.png +W1000 +H1000 +A");

