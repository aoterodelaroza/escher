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

function a = mol_angle(mol, at1, at2, at3)
% function a = mol_angle(mol, at1, at2, at3)
%
% mol_angle - returns the angle (in degrees) between three atoms in the molecule. 
%
% Required input variables:
% mol: structure containing the molecule.
% at1, at2, at3: return the at1-at2-at3 angle.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: Nov 2012

  x21 = mol.atxyz(:,at1) - mol.atxyz(:,at2);
  x23 = mol.atxyz(:,at3) - mol.atxyz(:,at2);
  a = acos(dot(x21,x23) / norm(x21) / norm(x23)) * 180 / pi;

endfunction
