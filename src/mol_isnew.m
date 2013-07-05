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

function isnew = mol_isnew(mol, xyz, eps=1e-8, LOG=1)
% function isnew = mol_isnew(mol, xyz, eps=1e-8, LOG=1)
%
% mol_isnew - tests if the atom of coordinates xyz is closer than eps to
%    any atom in mol.
%
% Required input variables:
% {mol}: structure containing the molecule.
% {xyz}: vector with the coordinates to test.
%
% Optional input variables (all have default values):
% {eps}: distance to consider the atoms degenerate.
% {LOG}: print detailed information if LOG>0.
%
% Output variables:
% {isnew}: result of the test.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2012

j = 1;
isnew = true;
while (j<=length(mol.atnumber) && isnew)
   if (size(xyz)==[3,1])
      x = xyz - mol.atxyz(1:3,j);
   elseif (size(xyz)==[3,1])
      x = xyz' - mol.atxyz(1:3,j);
   else
      error('mol_isuniq: wrong xyz dimensions!');
   endif
   d = sqrt(x' * x);
   if (d <= eps)
      if (LOG>0)
         printf('mol_isnew: coordinates correspond to atom %d in mol\n', j);
      endif
      isnew = false;
      return
   endif
   j++;
endwhile
isnew = true;
if (LOG>0)
   printf('mol_isnew: The atom is new!\n');
endif
endfunction
