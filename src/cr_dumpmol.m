% Copyright (C) 2015 Alberto Otero-de-la-Roza and Victor Lua~na
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

function mol = cr_dumpmol(cr);
% function mol = cr_dumpmol(cr);
%
% cr_dumpmol - Write all the atoms in the unit cell to a molecule, 
% in Cartesian coordinates that are consistent with the Cartesian 
% representation of the lattice vectors. This routine is used to 
% modify the cell motif, along with cr_readmol().
%
% Input:
% cr: input crystal.
%
% Output:
% mol: a molecule containing all the atoms in the unit cell.
% The Cartesian coordinates are consistent with the crystal's
% lattice vectors.
%

  mol = molecule_();
  for i = 1:cr.nat
    mol = mol_addatom(mol,cr.attyp{cr.typ(i)},cr_x2c(cr,cr.x(i,:))*0.52917720859);
  endfor

endfunction
