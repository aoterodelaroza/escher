#! /etc/alternatives/octave -q

list = {"c6h6_c6h6_stack","ch4_ch4","ch4_hf",...
        "h2co_h2co","h2o_h2o","phenol_phenol"};
rlist = cell(1,length(list));
rlist{1} = rlist{4} = rlist{6} = ...
    [
     1 0 0   
     0 1 0   
     0 0 1   
     0 0 -15 
     ];

for i = 1:length(list)
  mol = mol_readxyz(sprintf("%s.xyz",list{i}));
  rep = mol_ball(mol)
  rep = mol_stick(mol,rep);
  rep = mol_stick(mol,rep,"O","H",[1.8,2.0],:,0.01,[0 255 0]);
  if (!isempty(rlist{i}))
    rep = rep_setdefaultscene(rep,rlist{i});
  else
    rep = rep_setdefaultscene(rep);
  endif
  rep_write_pov(rep,sprintf("%s.pov",list{i}));
  system(sprintf("povray -D -UV +I%s.pov +O%s.png +W1000 +H1000 +A",list{i},list{i}));
endfor
