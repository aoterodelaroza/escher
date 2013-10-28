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

function [cm,mass] = mol_cmass(mol, iat=[], LOG=1)
% function [cm,mass] = mol_cmass(mol, iat=[], LOG=1)
%
% mol_cmass - determine the center of mass of a molecule or part.
%
% Required input variables:
% {mol}: the molecule.
%
% Optional input variables (all have default values):
% {iat}: list of atoms in the molecule that will be taken into account.
%     By default (an empty iat) use the whole molecule.
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {cm}: column vector with the center of mass coordinates.
% {mass}: total mass of the molecule or fragment.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

# Check input data, particularly iat1 and iat2.
if (isempty(iat))
   iat = 1:n;
endif

# Get the center of mass:
xi = mol.atxyz(1:3,iat); massi = mol.atmass(iat);
mass = sum(massi);
cm = xi * mol.atmass(iat)' / mass;

if (LOG > 0)
   printf('mol_cmass: center of mass of (part of?) a molecule\n');
   printf('Fragment (%d):', length(mol.atmass));
   printf(' %d', iat);
   printf('\n');
   for i = 1 : length(iat)
      ii = iat(i);
      printf('%4d %2s', ii, mol.atname{ii});
      printf('%12.6f%12.6f%12.6f', mol.atxyz(1:3,ii));
      printf('\n');
   endfor
   printf('\nCenter of mass:');
   printf('%12.6f%12.6f%12.6f', cm);
   printf('\n');
endif

endfunction
