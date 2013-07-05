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

function molout = mol_addatom (atname, atxyz, molin, newmol=0, LOG=1)
% function molout = mol_addatom (atname, atxyz, molin, newmol=0, LOG=1)
%
% mol_addatom - add a new atom (if it is not included) to the molecule.
%
% Required input variables:
% atmane: name of new atom.
% atxyz: (1:3) array with the cartesian coordinates of the new atom.
% molin: structure with the input molecular description. The format is:
%       molin.name ---> name of the molecule.
%       molin.atname --> {1:M} cell array with the symbols of the atoms
%                        (M is the number of atoms in the molecule).
%       molin.atxyz ----> Mx3 matrix with the atomic coordinates.
%       molin.atmass --> [1:M] vectos with atomic masses.
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
% {newmol = 0}: Enter newmol!= to create a new molecule or clean and restart
%               an old molecule.
%
% Required output variables:
% molout: structure with the input molecular description. The format is:
%       molout.name ---> name of the molecule.
%       molout.atname --> {1:M} cell array with the symbols of the atoms
%                        (M is the number of atoms in the molecule).
%       molout.atxyz ---> Mx3 matrix with the atomic coordinates.
%       molout.atmass --> [1:M] vectos with atomic masses.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

  angtobohr = 1.88972613288564;
  bohrtoans = 0.52917720859;

  if (newmol != 0)
     mol.name = 'unknown';
     mol.atname = {};
     mol.atnumber = [];
     mol.atxyz = [];
     mol.atmass = [];
     new = 1;
  else
     mol = molin;
     new = length(mol.atname) + 1;
  endif

  mol.atname{new} = atname;
  [atnumber,atprop] = mol_dbatom(atname);
  mol.atnumber(new) = atnumber;
  mol.atxyz(1:3,new) = atxyz(1:3);
  mol.atmass(new) = atprop.mass;

  molout = mol;

endfunction
