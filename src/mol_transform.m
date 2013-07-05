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

function molout = mol_transform (op, molin, LOG=1)
% function molout = mol_transform (op, molin, LOG=1)
%
% mol_transform - apply the "op" 3x4 matrix to the coordinates of "molin".
%
% Required input variables:
% op: 3x4 matrix containing the operation.
% molin: structure with the input molecular description. The format is:
%       molin.name ---> name of the molecule.
%       molin.atname --> {1:M} cell array with the symbols of the atoms
%                        (M is the number of atoms in the molecule).
%       molin.xyz -----> Mx3 matrix with the atomic coordinates.
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
%
% Required output variables:
% molout: structure with the input molecular description. The format is:
%       molout.name ---> name of the molecule.
%       molout.atname --> {1:M} cell array with the symbols of the atoms
%                        (M is the number of atoms in the molecule).
%       molout.xyz -----> Mx3 matrix with the atomic coordinates.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

  mol = molin;
  mol.atxyz = op(1:3,1:3) * mol.atxyz + op(1:3,4);
  molout = mol;

endfunction
