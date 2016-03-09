% Copyright (C) 2011-2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [ivect,ival] = mol_inertiamatrix (mol, LOG=0)
% function [ivect,ival] = mol_inertiamatrix (mol, LOG=0)
%
% mol_inertiamatrix - calculate and diagonalize the inertia matrix.
%
% Required input variables:
% {mol}: descrition of the molecule containing, at least:
%    mol.name ---> name of the molecule.
%    mol.atname --> {1:M} cell array with the symbols of the atoms
%                   (M is the number of atoms in the molecule).
%    mol.atnumber-> [1:M] vector of atomic numbers.
%    mol.atxyz ---> Mx3 matrix with the atomic coordinates.
%    mol.atmass --> [1:M] vector with atomic masses.
%
% Optional input variables (all have default values):
% {LOG=0}: print information about the data read in if LOG>0.
%
% Required output variables:
% {ivect,ival}: eigenvectors (columns) and eigenvalues of the inertia matrix.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: January 2012

imat = zeros(3,3);
for k = 1:mol.nat
   rr = mol.atxyz(1:3,k) * mol.atxyz(1:3,k)';
   for i = 1:3
      imat(i,i) += mol.atmass(k) * (trace(rr) - rr(i,i));
      for j = 1:i-1
         imat(i,j) -= mol.atmass(k) * rr(i,j);
         imat(j,i) -= mol.atmass(k) * rr(j,i);
      endfor
   endfor
endfor
[ivect,ival] = eig(imat);
ival = diag(ival);

if (LOG>0)
   printf('mol_inertiamatrix (%s):\n', mol.name);
   printf('... Inertia matrix:\n');
   printf('    %15.9f %15.9f %15.9f\n', imat);
   printf('... Eigenvalues:\n');
   printf('    %15.9f %15.9f %15.9f\n', ival);
   printf('... Eigenvectors (columns):\n');
   printf('    %15.9f %15.9f %15.9f\n', ivect);
endif

endfunction
