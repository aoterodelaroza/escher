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

function mol = mol_getfragment (moli, iat, LOG=0)
% function mol = mol_getfragment (moli, iat, LOG=0)
%
% mol_getfragment - create a fragment from a list of atoms in a molecule.
%
% Required input variables:
% {moli}: input molecule
% {iat}: vector with the atoms to be copied from mol to the fragment.
%
% Optional input variables (all have default values):
% {LOG=0}: print information about the data read in if LOG>0.
%
% Required output variables:
% {mol}: fragment description
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2012

if (length(iat) < 1)
   error('mol_getfragment: iat list empty!');
elseif (min(iat)<1 || max(iat)>moli.nat)
   error('mol_getfragment: wrong iat list!');
endif

nm = sprintf(' %d', iat);
if (isfield(moli,"name") && !isempty(moli.name))
   mol.name = sprintf('Fragment[%s] of %s', nm, moli.name);
endif
mol.nat = length(iat);
mol.atname = moli.atname(iat);
mol.atnumber = moli.atnumber(iat);
mol.atxyz = moli.atxyz(:,iat);
mol.atmass = moli.atmass(iat);

if (LOG>0)
   printf('moli_getfragment:\n');
   if (exist(mol.name))                             
     printf('... %s\n', mol.name);
   endif
   printf("number, orig, symbol, at_number, xyz coordinates:\n");
   for i = 1:length(iat)
      printf("%3d %3d %2s %3d", i, iat(i), mol.atname{i}, mol.atnumber(i));
      printf(" %12.6f %12.6f %12.6f", mol.atxyz(1:3,i));
      printf("\n");
   endfor
endif

endfunction
