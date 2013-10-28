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

function [check,vector] = mol_islinear(mol, eps=1e-6, LOG=0)
%function [check,vector] = mol_islinear(mol, eps=1e-6, LOG=0)
%
% mol_islinear - check if the molecule is linear (within a given eps
%    threshold) and, if it is, determines the vector line.
%
% Required input variables:
% {mol}: the molecule.
%
% Optional input variables (all have default values):
% {eps}: accepted error on the linearity test.
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {check}: true/false linearity result.
%
% Optional output variable:
% {vector}: normalized vector describing the molecule orientation from the
%           center of mass.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2012

# Center of mass:
n = mol.nat;
cm = mol.atxyz * mol.atmass' / sum(mol.atmass);
# Place cm at the origin:
x = mol.atxyz - cm * ones(1,n);

# The fourth component, x(4,i), will be the modulus of x(1:3,i) (i = 1:n)
# Be careful with atoms at the cm:
for i = 1:n
   x(4,i) = sqrt( dot( x(1:3,i), x(1:3,i) ) );
endfor

# Test for linear molecules:
# The normalized scalar product of the position vectors of each pair
# of atoms must be +1 or -1.
# (adapted from tessel2).
summ = 0;
for i = 1:n-1
   if (x(4,i) > eps)
      for j = i+1:n
         if (x(4,j) > eps)
            pesc = (x(1:3,i)' * x(1:3,j)) / (x(4,i) * x(4,j));
            contri = abs(pesc) - 1;
            summ += abs(contri);
            if (summ > eps)
               check = false; vector = [0;0;0];
               if (LOG > 0)
                  printf('mol_islinear: This molecule is NOT linear!\n');
               endif
               return
            endif
         endif
      endfor
   endif
endfor

# Molecule axis vector (normalized):
# Use the atom most distant from cm to reduce errors.
[xmax,imax] = max(x(4,:));
if (xmax > eps)
   vector = x(1:3,imax) / xmax;
else
   error('mol_linear: atoms are too close to cm!');
endif
check = true;

if (LOG > 0)
   printf('mol_islinear: test result -> %d\n', check);
   printf('... linearity sum (must be below %.2e) -> %.3e\n', eps, summ);
   if (check)
      printf('... center of mass --> %.6f %.6f %.6f\n', cm);
      printf('... mol axis vector -> %.6f %.6f %.6f\n', vector);
   endif
endif

endfunction
