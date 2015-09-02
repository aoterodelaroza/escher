% Copyright (c) 2015 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [form, mass] = mol_formula (mol, LOG=1)
% function [form, mass] = mol_formula (mol, LOG=1)
%
% mol_formula - determine the molecular formula of mol
%
% Required input variables:
% mol: cell array with the molecular structure data.
%      This is automatically created by mol_readcube.m and similar
%       routines.
%
% Essential elements in the mol cell array:
% mol.nat ... number of atoms
% mol.atname .... cell array with the atom names
% mol.atnumber .. array with the integer atomic numbers
%
  [at_types at_typ_num j] = unique(mol.atname);
  g = '';
  for k = length(at_types):-1:1
    g = sprintf("%s%s%d",g,mol.atname{at_typ_num(k)},sum(j == k));
  endfor
  form = g;

  mass = 0;
  for k = 1:length(at_types)
    symb = at_types(k);
    [Z,at] = mol_dbatom(symb);
    mass += at_typ_num(k) * at.mass;
  endfor

endfunction
