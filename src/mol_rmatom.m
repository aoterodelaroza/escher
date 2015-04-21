% Copyright (C) 2011-2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function mol = mol_rmatom(molin, id)
% function mol = mol_rmatom(molin, id)
%
% mol_rmatom - remove a list of atoms from the molecule.
%
% Input:
% molin: input molecule.
% id: integer list of atoms.
%
% Required output variables:
% mol: output molecule.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

  mol = molin;
  ikeep = setdiff(1:mol.nat,id);
  mol.nat = length(ikeep);
  mol.atname = mol.atname(ikeep);
  mol.atnumber = mol.atnumber(ikeep);
  mol.atmass = mol.atmass(ikeep);
  mol.atxyz = mol.atxyz(:,ikeep);

endfunction
