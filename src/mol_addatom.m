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

function mol = mol_addatom (molin, atname, atxyz, LOG=0)
% function mol = mol_addatom (molin, atname, atxyz, LOG=0)
%
% mol_addatom - add a new atom (if it is not included) to the molecule.
%
% Required input variables:
% atname: name of new atom.
% atxyz: (1:3) array with the cartesian coordinates of the new atom.
% molin: structure with the input molecular description. 
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
%
% Required output variables:
% mol: structure with the input molecular description. The format is:
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: June 2011

  angtobohr = 1.88972613288564;
  bohrtoans = 0.52917720859;

  mol = molin;
  mol.nat = mol.nat + 1;

  mol.atname{mol.nat} = atname;
  [atnumber,atprop] = mol_dbatom(atname);
  mol.atnumber(mol.nat) = atnumber;
  mol.atxyz(1:3,mol.nat) = atxyz(1:3);
  mol.atmass(mol.nat) = atprop.mass;

endfunction
