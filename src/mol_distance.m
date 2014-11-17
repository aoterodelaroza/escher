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

function d = mol_distance (mol, i1, i2)
% function d = mol_distance (mol, i1, i2)
%
% mol_distance - calculate the distance (in angstrom) between two
% atoms (i1 and i2) in the molecule.
%
% Required input variables:
% mol: the input molecule
% i1, i2: index for the atoms

  d = norm(mol.atxyz(:,i1)-mol.atxyz(:,i2));

endfunction
