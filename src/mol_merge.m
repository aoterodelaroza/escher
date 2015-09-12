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

function [mol] = mol_merge(varargin);
% function [mol] = mol_merge(mol1, mol2,...)
%
% mol_merge - merge molecules.
%
% Required input variables:
% {mol1,mol2,...}: molecules or cell arrays of molecules to merge.
%
% Required output variables:
% {mol}: merged molecule. 
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

mol = molecule_();
mol.name = "merged: ";
n = 0;
for i = 1:length(varargin)
  if (!iscell(varargin{i}))
    mcell = {varargin{i}};
  else
    mcell = varargin{i};
  endif

  for k = 1:length(mcell)
    mol1 = mcell{k};
    n1 = mol1.nat;
    if (isfield(mol1,"name") && !isempty(mol1.name))
      mol.name = sprintf("%s %s",mol.name,mol1.name);
    endif
    for j = 1 : n1
      n++;
      mol.atxyz(1:3,n) = mol1.atxyz(1:3,j);
      mol.atname{n} = mol1.atname{j};
      mol.atnumber(n) = mol1.atnumber(j);
      mol.atmass(n) = mol1.atmass(j);
    endfor
  endfor
endfor
mol.nat = n;

endfunction
