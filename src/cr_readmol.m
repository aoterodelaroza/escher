% Copyright (C) 2015 Alberto Otero-de-la-Roza and Victor Lua~na
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

function cr = cr_readmol(cr0,mol);
% function cr = cr_readmol(cr0,mol);
%
% cr_readmol - Read the atoms from the molecule mol into the crystal,
% replacing the current molecular motif. The crystal cell is preserved
% and the input lattice vectors are used to convert the Cartesian
% coordinates in the molecule to crystallographic coordinates. This routine
% together with cr_dumpmol is useful for manipulating the crystal
% motif with an external molecular editor.
%
% Input:
% cr: input crystal.
% mol: input molecule. The molecule is read and replaces all the atoms
% in the molecular motif of cr. The lattice vectors are preserved.
%
% Output:
% mol: a crystal with the same lattice vectors as the input but with the 
% molecular motif contained in mol.
%

  cr = crystal_();
  cr.name = cr0.name;
  cr.a = cr0.a;
  cr.b = cr0.b;
  cr.r = cr0.r;
  cr.g = cr0.g;
  cr.omega = cr0.omega;
  
  for i = 1:mol.nat
    cr = cr_addatom(cr,cr_c2x(cr,mol.atxyz(:,i)/0.52917720859),mol.atname{i});
  endfor

endfunction
