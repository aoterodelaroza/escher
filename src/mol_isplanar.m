% Copyright (C) 2011--12 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [check,vector] = mol_isplanar(mol, eps=1e-6, LOG=0)
%function [check,vector] = mol_isplanar(mol, eps=1e-6, LOG=0)
%
% mol_isplanar - check if the molecule is planar (within a given eps
%    threshold) and, if it is, determine the plane normal vector.
%
% Required input variables:
% {mol}: the molecule.
%
% Optional input variables (all have default values):
% {eps}: accepted error on the planarity test.
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {check}: true/false planarity result.
%
% Optional output variable:
% {vector}: normalized vector describing the normal to the molecular plane.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: January 2012

# Center of mass:
n = mol.nat;
cm = mol.atxyz * mol.atmass' / sum(mol.atmass);
# Place cm at the origin:
x = mol.atxyz - cm * ones(1,n);

# Test for planar molecules:
# The triple (volume) product of the position vectors of each group
# of three different atoms must be zero.
# (adapted from tessel2).
summ = 0;
for i = 1:n-2
   for j = i+1:n-1
      for k = j+1:n
         pmixto = dot(x(1:3,i), cross(x(1:3,j), x(1:3,k)));
         summ += abs(pmixto);
         if (summ > eps)
            check = false; vector = [0;0;0];
            if (LOG > 0)
               printf('mol_isplanar: This molecule is NOT planar!\n');
               printf('... offending trio %d %d %d (%.3e)\n', i, j, k, summ);
            endif
            return
         endif
      endfor
   endfor
endfor

# Plane normal vector (normalized):
# Use the atom couple with the largest cross vector to reduce errors.
pv = zeros(4,(n*(n-1))/2);
m = 0;
for i = 1:n-1
   for j = i+1:n
      pv(1:3,++m) = cross(x(1:3,i), x(1:3,j));
      pv(4,m) = dot(pv(1:3,m), pv(1:3,m));
   endfor
endfor
[xmax,imax] = max(pv(4,:));
if (xmax > eps)
   vector = pv(1:3,imax) / sqrt(xmax);
else
   error('mol_planar: atoms are too close to cm!');
endif
check = true;

if (LOG > 0)
   printf('mol_isplanar: test result -> %d\n', check);
   printf('... planarity sum (must be below %.2e) -> %.3e\n', eps, summ);
   if (check)
      printf('... center of mass --> %.6f %.6f %.6f\n', cm);
      printf('... normal vector ---> %.6f %.6f %.6f\n', vector);
   endif
endif

endfunction
