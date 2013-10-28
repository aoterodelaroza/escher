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

function D = mol_dist2(mol1, iat1, mol2, iat2, LOG=0)
% function D = mol_dist2(mol1, iat1, mol2, iat2, LOG)
%
% mol_dist2 - returns the matrix of distances between the list iat1 of atoms
%     in the molecule mol1, and the list iat2 of atoms in the molecule mol2.
%
% Required input variables:
% {mol1,mol2}: structures containing the molecules.
% {iat1,iat2}: lists of atoms whose distance is being calculated.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2012

D = zeros(length(iat1),length(iat2));
for i = 1 : length(iat1)
   ii = iat1(i);
   for j = 1 : length(iat2)
      jj = iat2(j);
      x = mol2.atxyz(1:3,jj) - mol1.atxyz(1:3,ii);
      D(i,j) = sqrt(x' * x);
   endfor
endfor

if (LOG > 0)
   printf('mol_dist2:\n');
   printf('Atoms in mol1 ->'); printf(' %d', iat1); printf('\n');
   printf('Atoms in mol2 ->'); printf(' %d', iat2); printf('\n');
   printf('Matrix of distances:\n');
   for i = 1 : length(iat1)
      for j = 1 : length(iat2)
         printf('%12.5f', D(i,j));
      endfor
      printf('\n');
   endfor
endif

endfunction
