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

function geom = mol_groupgeometry (mol, g1, g2, LOG=1)
% function geom = mol_groupgeometry (mol, g1, g2, LOG=1)
%
% mol_groupgeometry - Given a molecular database, mol, and two groups made
% of subsets of mol, this routine calculates and prints several geometry
% parameters of g1 with respect to g2.
%
% Required input variables:
% mol: molecular database, containing, at least:
%      mol.atxyz[1:3,1:nat] ...array of cartesian coordinates.
% g1[:],g2[:]: index of the atoms in mol that belong to g1 and g2.
%
% % Optional input variables (all have default values):
% {LOG}: current level of printing. 0: no printing.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: November 2011

% 1st test: distances between g1 and g2:
dist = zeros(length(g1),length(g2));
for i = 1 : length(g1)
   for j = 1 : lenght(g2)
      x = mol.atxyz(1:3,g1(i)) - mol.atxyz(1:3,g2(j));
      dist(i,j) = sqrt(dot(x,x));
   endfor
   dist1 = dist(i,:);
   [dm,im] = min(dist1);
   printf('Closest to g1(%d) in g2: (%d) --> %.6f\n', \
         g1(i), g2(im), dm);
endfor
printf('\n');
for j = 1:length(g2)
   dist1 = dist(:,j);
   [dm,im] = min(dist1);
   printf('Closest to g2(%d) in g1: (%d) --> %.6f\n', \
         g2(j), g1(im), dm);
endfor


endfunction
