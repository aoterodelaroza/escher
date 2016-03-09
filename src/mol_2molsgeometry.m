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

function void = mol_2molsgeometry (mol1, mol2, LOG=1)
% function void = mol_2molsgeometry (mol1, mol2, LOG=1)
%
% mol_2molsgeometry - Given a couple of molecular databases, mol1 and
% mol2, this routine calculates and prints the distances between the atoms
% the first to the atoms in the second molecule.
%
% Required input variables:
% mol1: molecular database, containing, at least:
%       mol1.atname{1:nat} ... cell array of strings with the atom names.
%       mol1.atxyz[1:3,1:nat] ...array of cartesian coordinates.
% mol2: the same for the second molecule.
%
% % Optional input variables (all have default values):
% {LOG}: current level of printing. 0: no printing.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: October 2012

if (LOG>0)
   n1 = mol1.nat;
   n2 = mol2.nat;
   printf('mol_2molsgeometry:\n');
   printf('Atom(mol1)  Atom(mol2) dist Coords.\n');
   for i = 1 : n1
      dist = zeros(1,n2);
      for j = 1 : n2
         x = mol1.atxyz(1:3,i) - mol2.atxyz(1:3,j);
         dist(j) = sqrt(dot(x,x));
      endfor
      [dm,jm] = min(dist);
      printf('%4d %-5s ', i, cell2mat(mol1.atname(i)));
      printf('%4d %-5s %10.4f %10.4f %10.4f %10.4f', ...
             jm, cell2mat(mol2.atname(jm)), dm, mol2.atxyz(1:3,jm));
      printf('\n');
      [ds,js] = sort(dist,2,'ascend');
      if (LOG > 1)
         for j = 2:min([LOG,n2,4])
            jm = js(j); dm = ds(j);
            printf('           %4d %-5s %10.4f %10.4f %10.4f %10.4f', ...
                   jm, cell2mat(mol2.atname(jm)), dm, mol2.atxyz(1:3,jm));
            printf('\n');
         endfor
         printf('\n');
      endif
   endfor
endif

endfunction
