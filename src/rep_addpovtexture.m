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

function rep_addpovtexture(name,string);
% function rep_addpovtexture(name,string)
%
% rep_addpovtexture - add a povray texture to the database.
%
% Required input variables:
% name: name of the new texture.
% string: povray texture definition.

  global texdb
  
  if (!exist("texdb","var") || isempty(texdb))
    rep_texdbstart();
  endif
  
  n = length(texdb);
  texdb{n+1}.typ = "pov";
  texdb{n+1}.name = name;
  texdb{n+1}.string = string;

endfunction
