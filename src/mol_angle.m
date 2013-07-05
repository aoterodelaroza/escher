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

function a = mol_angle(mol, at1, at2, at3, LOG=0)
% function a = mol_angle(mol, at1, at2, at3, LOG=0)
%
% mol_angle - returns the angle between three atoms in the molecule.
%
% Required input variables:
% mol: structure containing the molecule.
% at1, at2, at3: return the at1-at2-at3 angle.
%
% % Optional input variables (all have default values):
% {LOG}: current level of printing. 0: no printing.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Nov 2012

x21 = mol.atxyz(1:3,at2) - mol.atxyz(1:3,at1);
x23 = mol.atxyz(1:3,at2) - mol.atxyz(1:3,at3);
pesc = dot(x21,x23);
pvec = cross(x21,x23);
pvec = sqrt(dot(pvec,pvec));
a = atan2(pvec,pesc)*180/pi;

if (LOG > 0)
   s1 = sprintf("%s(%d)", cell2mat(mol.atname(at1)), at1);
   s2 = sprintf("%s(%d)", cell2mat(mol.atname(at2)), at2);
   s3 = sprintf("%s(%d)", cell2mat(mol.atname(at3)), at3);
   printf("%12.6f %s-%s-%s\n", a, s1, s2, s3);
endif

endfunction
