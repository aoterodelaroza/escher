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

function molout = mol_transform (op, t, molin, LOG=1)
% function molout = mol_transform (op, t, molin, LOG=1)
%
% mol_transform - apply the "op" rotation (3x3) and the t translation (3x1) to the coordinates of "molin".
%
% Required input variables:
% op: 3x3 matrix containing the rotation.
% t: 3x1 translation vector
% molin: structure with the input molecular description. 
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
%
% Required output variables:
% molout: structure with the input molecular description. 
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

  mol = molin;
  mol.atxyz = op * mol.atxyz + t;
  molout = mol;

endfunction
