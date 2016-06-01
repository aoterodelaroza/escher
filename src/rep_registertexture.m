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

function [rep itex] = rep_registertexture(rep0,tex)
% function [rep itex] = rep_registertexture(rep,tex)
%
% rep_registertexture - register the texture tex (a string, corresponding
% to an entry in the texdb) in the representation. 
%
% Input:
% rep: the input representation.
% tex: string for the texture.
%
% Output:
% rep: the representation containing the new texture.
% itex: integer index for the new texture.
%

  global texdb

  ## See if we know this texture already
  rep = rep0;
  anames = fieldnames(rep.texlib);
  for i = 1:length(anames)
    name = anames{i};
    if (strcmpi(tex,name))
      itex = name;
      return
    endif
  endfor

  ## Initialize the texture library, if not initialized already
  if (!exist("texdb","var") || isempty(texdb))
    tex_dbstart();
  endif

  ## Search for the texture in the database
  found = 0;
  for i = 1:length(texdb)
    if (strcmpi(tex,texdb{i}.name))
      found = i;
      break
    endif
  endfor

  if (found > 0)
    itex = tex;
    rep.texlib = setfield(rep.texlib,tex,found);
  else
    error("texture not found")
  endif

endfunction
