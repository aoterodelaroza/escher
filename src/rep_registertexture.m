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

  rep = rep0;
  n = length(rep.texlib);
  for i = 1:n
    if (strcmp(tex,rep.texlib{i}))
      itex = i;
      return
    endif
  endfor
  itex = n+1;
  rep.texlib{n+1} = tex;

endfunction
