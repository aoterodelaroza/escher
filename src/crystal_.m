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

function cr = crystal_()
% function cr = crystal_()
%
% crystal - create an empty crystal structure and initialize the number
% of atoms to zero. Cell array constructor for an empty crystal.
%
% Output:
% {cr}: the empty crystal structure with all the fields defined.
%

  cr.name = "";
  cr.nat = 0;               # Number of atoms in the unit cell
  cr.ntyp = 0;              # Number of atom types
  cr.attyp = cell();        # Atomic symbols
  cr.rvdwtyp = cr.c6typ = cr.zvaltyp = cr.qtyp = cr.ztyp = [];
     # Parameters for each atomic type
     # cr.rvdwtyp .... van der Waals radius (bohr)
     # cr.c6typ ...... virial C6 coeficient
     # cr.zvaltyp .... number of valence electrons
     # cr.qtyp ....... Atomic charge
     # cr.ztyp ....... Atomic number (Z)
  cr.typ = [];              # Atom type indexes (1..ntyp). The symbol of typ 1
                            # would reside in cr.attyp(cr.typ(1))
  cr.x = [];                # Atom crystal coordinates
  cr.a = cr.b = zeros(1,3); # unit cell lengths (a,bohr) and angles (b,radians)
                            # That's the basic definition of the crytal parallelepiped
  cr.g = zeros(3);          # Real space metric matrix
  cr.r = zeros(3);          # Crystallographic to cartesian matrix (bohr)
  cr.omega = 0;             # Unit cell volume

endfunction
