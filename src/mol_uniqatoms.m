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

function [uniqlist,excludelist] = mol_uniqatoms(mol, eps=1e-8, LOG=1)
% function [uniqlist,excludelist] = mol_uniqatoms(mol, eps=1e-8, LOG=1)
%
% mol_uniqatoms - detects atoms closer than eps and returns a list of the
%    unique atoms.
%    Notice: only the distances are considered!!!
%
% Required input variables:
% {mol}: structure containing the molecule.
%
% Optional input variables (all have default values):
% {eps}: distance to consider the atoms degenerate.
% {LOG}: print detailed information if LOG>0.
%
% Output variables:
% {uniqlist}: list of unique atoms.
%
% Optional output variables:
% {excludelist}: list of atoms that could be a copy of the unique list.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2012

nuniq = 0;
nexclude = 0;
excludelist = [];
uniqlist(++nuniq) = 1;
n = length(mol.atnumber);
for i = 2:n
   isuniq = true;
   j = 1;
   while (j<=nuniq && isuniq)
      x = mol.atxyz(1:3,i) - mol.atxyz(1:3,j++);
      d = sqrt(x' * x);
      if (d <= eps)
         isuniq = false;
      endif
   endwhile
   if (isuniq)
      uniqlist(++nuniq) = i;
   else
      excludelist(++nexclude) = i;
   endif
endfor

if (LOG>0)
   printf('mol_uniqatoms:\n',nuniq,nexclude);
   printf('Unique (%d) and repeated (%d) atoms\n', nuniq, nexclude);
   printf('Unique atoms --->'); printf(' %d', uniqlist);    printf('\n');
   if (nexclude > 0)
      printf('Repeated atoms ->'); printf(' %d', excludelist); printf('\n');
   else
      printf('Repeated atoms -> empty list!\n');
   endif
endif

endfunction
