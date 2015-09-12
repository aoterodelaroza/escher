% Copyright (c) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
% Adapted from a tessel routine.
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

function [mol] = cr_applymask(cr, mask);
% function [mol] = cr_applymask(cr, mask);
%
% cr_applymask - obtains a molecule from a crystal by applying a crystal mask.
%
% Required input variables:
% {cr}: struct containing the crystal.
% {mask}: the crystal mask to apply.
%
% Required output variables:
% {mol}: the molecular description.
%

bohr2angstrom = 0.52917720859;

r = cr.r;

mol = molecule_();
mol.atname = cell(1,mask.nat);
mol.atnumber = zeros(1,mask.nat);
mol.atxyz = zeros(3,mask.nat);
n = 0;
for i = 1:mask.nat
  k = mask.i(i);
  x = cr.x(k,:) + mask.l(i,:);
  n += 1;
  mol.atname{n} = cr.attyp{cr.typ(k)};
  mol.atnumber(n) = cr.ztyp(cr.typ(k));
  mol.atxyz(:,n) = x * r;
endfor

if (n > 0)
  mol.atxyz = mol.atxyz(:,1:n) * bohr2angstrom;
else
  mol.atnumber = mol.atxyz = [];
endif

if (isfield(cr,"name") && !isempty(cr.name))
  mol.name = cr.name;
endif
mol.nat = n;
mol = mol_fillatmass(mol);

endfunction
