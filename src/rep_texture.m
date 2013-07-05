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

function tex = rep_texture(typ,name);
% function tex = rep_texture(typ,name);
%
% rep_texture - return one of the textures in the internal texture
% database, for using in povray or obj files.
%
% Required input variables:
% typ: type of texture, "pov" or "povray" for povray and
%      "obj" for obj.
% name: internal database handle of the texture. It has to be either
%       one of the default textures or one added using rep_addtexture.

  global texdb
  
  if (!exist("texdb","var") || isempty(texdb))
    rep_texdbstart();
  endif
  
  if (strcmp(lower(typ),"povray") || strcmp(lower(typ),"pov"))
    typ = "pov";
  elseif (strcmp(lower(typ),"obj"))
    typ = "obj";
  else
    error("unrecognized texture format")
  endif

  ## search for the texture in the database
  found = 0;
  for i = 1:length(texdb)
    if (strcmp(texdb{i}.typ,typ) && strcmp(texdb{i}.name,lower(name)))
      tex = texdb{i}; 
      found = 1;
    endif
  endfor
  if (!found)
    printf("Warning: texture %s not found. Falling back to tdefault.\n",name);
    tex = rep_texture(typ,"tdefault");
  endif

endfunction
