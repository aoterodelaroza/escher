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

function molout = mol_unitconvert (mol, conv='bohr2angstrom', LOG=1)
% function molout = mol_unitconvert (mol, conv='bohr2angstrom', LOG=1)
%
% mol_unitconvert - conversion of units in the mol database.
%
% Required input variables:
% {mol}: structure with the input molecular description. The format is:
%
% Optional input variables (all have default values):
% {conv}: conversion to apply. Possible conversions:
%    * bohr2angstrom
%    * angstrom2bohr
% {LOG = 1}: print information about the data read in if LOG>0.
%
% Required output variables:
% molout: structure with the output molecular description.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2012

angtobohr = 1.88972613288564;
bohrtoans = 0.52917720859;

molout = mol;
switch conv
   case {"bohr2angstrom"}
      molout.atxyz = mol.atxyz .* bohrtoans;
   case {"angstrom2bohr"}
      molout.atxyz = mol.atxyz .* angtobohr;
   otherwise
      error('mol_unitconversion: invalid conversion!');
endswitch

if (LOG > 0)
   printf('mol_unitconversion: %s\n', conv);
   printf('Final result:\n');
   natoms = molout.nat;
   printf("Number of atoms: %d\n", natoms);
   printf("number, symbol, at_number, xyz coordinates:\n");
   for i = 1:natoms
      printf("%4d %2s %3d", i, molout.atname{i}, molout.atnumber(i));
      printf(" %12.6f %12.6f %12.6f", molout.atxyz(1:3,i));
      printf("\n");
   endfor
endif

endfunction
