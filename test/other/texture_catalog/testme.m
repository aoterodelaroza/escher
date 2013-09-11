#! /etc/alternatives/octave -q

addpath("../../../src");

global texdb
rep_texdbstart();

rep = representation();
for i = 1:length(texdb)
  if (strcmp(texdb{i}.typ,"pov"))
    n = rep.nball = rep.nball + 1;
    rep.ball{n}.x = [2*(i-1), 0, 0];
    rep.ball{n}.r = 0.9;
    rep.ball{n}.rgb = [255 0 0];
    rep.ball{n}.tex = texdb{i}.name;
    n = rep.nball = rep.nball + 1;
    rep.ball{n}.x = [2*(i-1), 3, 0];
    rep.ball{n}.r = 0.9;
    rep.ball{n}.rgb = [0 255 0];
    rep.ball{n}.tex = texdb{i}.name;
    n = rep.nball = rep.nball + 1;
    rep.ball{n}.x = [2*(i-1), 6, 0];
    rep.ball{n}.r = 0.9;
    rep.ball{n}.rgb = [0 0 255];
    rep.ball{n}.tex = texdb{i}.name;
  endif
endfor
r = [
     1 0 0 
     0 1 0 
     0 0 1 
     0 0 -30
     ];

# rep_write_obj(rep,"textures.obj");
rep = rep_setdefaultscene(rep,r);
rep_write_pov(rep,"textures.pov");
system("povray -D -UV +Itextures.pov +Otextures.png +W2000 +H2000 +A");
