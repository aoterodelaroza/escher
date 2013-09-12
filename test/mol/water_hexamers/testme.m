#! /usr/bin/octave -q

addpath("../../../src");

list={"bag","boat1","boat2","book1","book2","cage","chair","prism"};
# list={"cage"};

for i = 1:length(list)
  mol = mol_readxyz(sprintf("%s.xyz",list{i}),0);
  rep = mol_ball(mol,:,"H",:,:);
  rep = mol_ball(mol,rep,"O",:,:);
  rep = mol_stick(mol,rep,"O","H",[0.8 1.1]);
  rep = mol_stick(mol,rep,"O","H",[1.5 2.3],:,0.01,[0 0 255]);
  rep_write_obj(rep,sprintf("%s.obj",list{i}));
  if (strcmp(list{i},"boat1"))
    r = [
         0.1691417694 0.6418759823 -0.7479215264 
         -0.9138480425 -0.1820913553 -0.3629391193
         -0.3691517711 0.7448747158 0.5557776093 
         0.0000000000 0.0000000000 -10.0000000000
         ];
    rep = rep_setdefaultscene(rep,r,0);
  elseif(strcmp(list{i},"boat2"))
    r = [
         0.0454561710 0.4583125710 -0.8876299262 
         -0.9467343688 -0.2637956142 -0.1846903861
         -0.3187965751 0.8487439156 0.4219061136 
         0.0000000000 0.0000000000 -10.0000000000
         ];
    rep = rep_setdefaultscene(rep,r,0);
  elseif(strcmp(list{i},"book2"))
    r = [
         1.0000000000 0.0000000000 0.0000000000 
         0.0000000000 1.0000000000 0.0000000000 
         0.0000000000 0.0000000000 1.0000000000 
         0.0000000000 0.0000000000 -10.0000000000
         ];
    rep = rep_setdefaultscene(rep,r,0);
  elseif(strcmp(list{i},"cage"))
    r = [
         0.7945815921 0.5708985925 -0.2066781074
         0.3985797167 -0.2336850166 0.8868657947
         0.4580122828 -0.7870650887 -0.4132304192
         0.0000000000 0.0000000000 -10.0000000000
         ];
    rep = rep_setdefaultscene(rep,r,0);
  elseif(strcmp(list{i},"chair"))
    r = [
         1.0000000000 0.0000000000 0.0000000000 
         0.0000000000 1.0000000000 0.0000000000 
         0.0000000000 0.0000000000 1.0000000000 
         0.0000000000 0.0000000000 -10.000000000
         ];
    rep = rep_setdefaultscene(rep,r,0);
  else
    r = [
         1.0000000000 0.0000000000 0.0000000000 
         0.0000000000 1.0000000000 0.0000000000 
         0.0000000000 0.0000000000 1.0000000000 
         0.0000000000 0.0000000000 -10.000000000
         ];
    rep = rep_setdefaultscene(rep,:,0);
  endif
  rep_write_pov(rep,sprintf("%s.pov",list{i}));
  system(sprintf("povray -D -UV +I%s.pov +O%s.png +W1000 +H1000 +A",list{i},list{i}));
endfor
