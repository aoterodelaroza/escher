% Copyright (C) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function rcov = mol_rcov(z)
% function rcov = mol_rcov(z)
%
% mol_rcov - returns the covalent radius of an atom with atomic number z
%
% Input:
% z: atomic number.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Nov 2012

  global atdb

  if (!exist("atdb","var") || isempty(atdb))
    err = mol_dbstart(LOG);
    if (err != 0)
      error("mol_dbatom: the atomic database does not start right!");
    endif
  endif
  
  if (z < 1 || z > length(atdb.rcov))
    rcov = 0;
  else
    rcov = atdb.rcov(z);
  endif

endfunction
