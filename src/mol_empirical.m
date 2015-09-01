% Copyright (c) 2015 Victor Lua~na and Alberto Otero-de-la-Roza
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

function g = mol_empirical (mol, LOG=1)
% function g = mol_empirical (mol, LOG=1)
%
% mol_empirical - determine the empirical formula of mol
%
% Required input variables:
% mol: cell array with the molecular structure data.
%      This is automatically created by mol_readcube.m and similar
%       routines.
%
% Essential elements in the mol cell array:
% mol.nat ... number of atoms
% mol.atname .... cell array with the atom names
% mol.atnumber .. array with the integer atomic numbers
%
mol.attyname = unique(mol.atname);
nt = length(mol.attname);
mol.attypen = zeros(nt);
g = '';
for t = 1:nt
   for k = 1:length(mol.atname)
      if (mol.attype(t) == mol.atname(k))
         mol.attypen(t) += 1;
      endif
   endfor
   g = sprintf('%s%d', mol.attype(t) , mol.attypen(t));
endfor

endfunction
