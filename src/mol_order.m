% Copyright (C) 2011--12 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [mol] = mol_order(moli, otype, mapfun="")
# function [mol] = mol_order(moli, otype, mapfun="")
%
% mol_order - reorder the atoms list.
%
% Required input variables:
% {moli}: the input molecule.
% {otype}: order by: 'atnumber' (atomic number). More to be added.
% 
% Optional input variables:
% {mapfun}: a function to apply to the order array. For instance, @fliplr to
%           reverse.
%
% Required output variables:
% {mol}: the output molecule.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: January 2012

if (ischar(otype)) 
  if (strcmpi(otype,'atnumber'))
    [s,i] = sort(moli.atnumber);
  endif
else (isnumeric(otype))
  i = otype;
endif

if (isa(mapfun,"function_handle"))
  i = mapfun(i);
endif
mol = moli;
mol.atname = moli.atname(i);
mol.atnumber = moli.atnumber(i);
mol.atxyz = moli.atxyz(:,i);
mol.atmass = moli.atmass(i);

endfunction
