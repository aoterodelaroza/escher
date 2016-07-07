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

function dh = mol_dihedral(mol, at1, at2, at3, at4)
% function dh = mol_dihedral(mol, at1, at2, at3, at4)
%
% mol_dihedral - returns the dihedral angle (in degrees) between four
% atoms in the molecule.
%
% Required input variables:
% mol: structure containing the molecule.
% at1, at2, at3, at4: return the at1-at2-at3-at4 dihedral angle (in degrees).
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: Nov 2012

  b1 = mol.atxyz(:,at2)-mol.atxyz(:,at1); 
  b2 = mol.atxyz(:,at3)-mol.atxyz(:,at2); 
  b3 = mol.atxyz(:,at4)-mol.atxyz(:,at3); 
  b12 = cross(b1, b2);
  b23 = cross(b2, b3);
  dh = atan2(norm(b2) * dot(b1, b23), dot(b12, b23)) * (180/pi);

endfunction
