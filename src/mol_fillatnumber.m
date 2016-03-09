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

function [mol] = mol_fillatnumber(moli)
% function [mol] = mol_fillatnumber(mol)
%
% mol_fillatnumber - fill in the atomic numbers in mol using the internal database
%                    and moli.atname.
%
% Required input variables:
% moli: the input molecule.
%
% Output variables:
% mol: molecule with atnumber array.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: June 2011

  mol = moli;
  for i = 1:moli.nat
    [Z,atprop] = mol_dbatom(moli.atname{i});
    mol.atnumber(i) = Z;
  endfor

endfunction
