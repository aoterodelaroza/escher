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

function [molout] = mol_reorder(molin, iorder1, iorder2, LOG=0)
%
% mol_reorder - change the order of atoms from the list in iorder1 to the
%    list in iorder2. For instance, mol_reorder(mol,1:6, 6:-1:1) reverses
%    the order of the first sixth atoms.
%
% Required input variables:
% {molin}: the input molecule.
% {iorder1}: initial order of the atoms.
% {iorder2}: final order of the atoms.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {molout}: the output molecule.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: January 2012

if (length(iorder1) != length(iorder2))
   error('mol_reorder: input and output orders have different lengths!');
endif

molout = molin;
molout.atname(iorder2) = molin.atname(iorder1);
molout.atnumber(iorder2) = molin.atnumber(iorder1);
molout.atxyz(:,iorder2) = molin.atxyz(:,iorder1);
molout.atmass(iorder2) = molin.atmass(iorder1);

if (LOG > 0)
   printf('mol_reorder:\n');
   printf('Initial order-->'); printf(' %d', iorder1); printf('\n');
   printf('Final order---->'); printf(' %d', iorder2); printf('\n');
endif

endfunction
