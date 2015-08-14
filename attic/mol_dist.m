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

function d = mol_dist(mol, at1, at2)
% function d = mol_dist(mol, at1, at2)
%
% mol_dist - returns the distance between two atoms in the molecule.
%
% Required input variables:
% mol: structure containing the molecule.
% at1, at1: atoms whose distance is being calculated.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

x = mol.atxyz(1:3,at1) - mol.atxyz(1:3,at2);
d = sqrt(x' * x);

endfunction
