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

function cr = cr_rmatom(cr0,id);
% function cr = cr_rmatom(cr0, id)
%
% cr_rmatom - remove atoms (or atoms) from the unit cell.
%
% Required input variables:
% cr0: input crystal.
% id: atom identifier (integer) or vector of atoms to remove.
%
% Output variables:
% cr: output crystal structure.
%

  if (!isscalar(id))
    alist = [id];
  else
    alist = id;
  endif
  
  idat = setdiff(1:cr0.nat,alist);
  cr = cr0;
  cr.nat = length(idat);
  cr.x = cr.x(idat,:);
  cr.typ = cr.typ(idat);

endfunction
