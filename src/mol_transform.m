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

function molout = mol_transform (molin, op, t=[0 0 0]')
% function molout = mol_transform (molin, op, t=[0 0 0]')
%
% mol_transform - apply the "op" rotation (3x3) and the t translation (3x1) to the coordinates of "molin".
%
% Required input variables:
% molin: structure with the input molecular description. 
% op: 3x3 matrix containing the rotation.
% t: 3x1 translation vector
%
% Optional input variables (all have default values):
%
% Required output variables:
% molout: structure with the input molecular description. 
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

  mol = molin;
  if (size(t) == [1 3])
    t = t';
  endif
  mol.atxyz = op * mol.atxyz + t * ones(1,mol.nat);
  molout = mol;

endfunction
