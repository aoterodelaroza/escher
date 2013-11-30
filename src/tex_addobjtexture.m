% Copyright (C) 2011 Victor Lua~na and Alberto Otero-de-la-Roza
%
% This octave routine is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version. See <http://www.gnu.org/licenses/>.
%
% The routine distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.

function tex_addobjtexture(name,Ns=96.078,Ka=[0 0 0],Ks=[128 128 128],Ni=1,illum=2);
% function tex_addobjtexture(name,Ns=96.078,Ka=[0 0 0],Ks=[0.5 0.5 0.5],Ni=1,illum=2);
%
% tex_addobjtexture - add an obj texture to the database.
%
% Required input variables:
% name: name of the new texture.
% Ns: shininess of the material.
% Ka: ambient color (three rgb integers)
% Ks: specular color (three rgb integers)
% Ni: optical density (index of refraction)
% illum: illumination model.

  global texdb
  
  if (!exist("texdb","var") || isempty(texdb))
    tex_dbstart();
  endif
  
  n = length(texdb);
  texdb{n+1} = texture();
  texdb{n+1}.typ = "obj";
  texdb{n+1}.name = name;
  texdb{n+1}.Ns = Ns;
  texdb{n+1}.Ka = Ka;
  texdb{n+1}.Ks = Ks;
  texdb{n+1}.Ni = Ni;
  texdb{n+1}.illum = illum;

endfunction
