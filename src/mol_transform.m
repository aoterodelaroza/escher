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

function mol = mol_transform (mol0, op, t=[0 0 0], local=0)
% function mol = mol_transform (molin, op, t=[0 0 0], local=0)
%
% mol_transform - apply the "op" rotation (3x3) and the t translation
% (3x1) to the coordinates of the input molecule.
%
% Required input variables:
% mol0: input molecule
% op: 3x3 matrix containing the rotation.
% t: 3x1 translation vector
% local: if true, apply the rotation at the molecular center of mass
% and not at the global coordinate origin.
%
% Optional input variables (all have default values):
%
% Required output variables:
% mol: output molecule.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: June 2011

  mol = mol0;
  if (all(size(t) == [1 3]))
    t = t';
  endif
  if (local)
    cm = mol_cmass(mol,1);
    mol.atxyz = mol.atxyz - cm * ones(1,mol.nat);
    mol.atxyz = op * mol.atxyz;
    mol.atxyz = mol.atxyz + cm * ones(1,mol.nat);
  else
    mol.atxyz = op * mol.atxyz;
  endif
  mol.atxyz = mol.atxyz + t * ones(1,mol.nat);

endfunction
