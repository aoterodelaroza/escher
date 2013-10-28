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

function d = mol_distmatrix(mol)
% function d = mol_distmatrix(mol)
%
% mol_distmatrix - returns the distance matrix of a molecule.
%
% Required input variables:
% {mol}: input molecule
%
% Output:
% d: distance matrix.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2012

M = mol.nat;
d = zeros(M);
for i = 1:M
  d(i,:) = norm(mol.atxyz - mol.atxyz(:,i) * ones(1,M),2,"columns");
endfor

endfunction
