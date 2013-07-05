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

function [mol] = mol_fillatmass(moli)
% function [mol] = mol_fillatmass(mol)
%
% mol_fillatmass - fill in the atomic masses in mol using the internal database
%                  and the moli.atnumber array.
%
% Required input variables:
% moli: the input molecule.
%
% Output variables:
% mol: molecule with atmass array.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

  mol = moli;
  for i = 1:length(moli.atnumber)
    [symb,atprop] = mol_dbsymbol(moli.atnumber(i));
    mol.atmass(i) = atprop.mass;
  endfor

endfunction
