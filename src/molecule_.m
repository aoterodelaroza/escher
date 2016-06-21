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

function mol = molecule_()
% function mol = molecule_()
%
% molecule - Coonstructor. Create an empty molecule structure and initialize mol.nat.
%
% Output:
% {mol}: the empty molecule structure with all the fields defined.
%        This can be the content of a unit cell in a crystal, a portion of the cell,
%        an isolated molecule, or even a fragment exceeding a unit cell.
%

  mol.nat = 0;
  mol.name = "";
  mol.atname = cell(); # List of the atom symbols in the molecule
  mol.atnumber = [];   # List of the Z atomic numbers
  mol.atmass = [];     # Atomic masses (gr/mol)
  mol.atxyz = [];      # Cartesian coordinates
  mol.adjl = [];       # Adjacency list (atomic connectivity)

endfunction
