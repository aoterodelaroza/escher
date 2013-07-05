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

function mol = cr_crystalbox(cr, x0=[-0.05 -0.05 -0.05], x1=[1.05 1.05 1.05], nmol=-1, LOG=0)
% function mol = cr_crystalbox(cr, x0=[-0.05 -0.05 -0.05], x1=[1.05 1.05 1.05], nmol=-1, LOG=0)
%
% cr_crystalbox - create a molecule by cutting a parallelepiped from a crystal.
%
% Required input variables:
% {cr}: struct containing the crystal.
%
% Optional input variables (all have default values):
% {x0}: one vertex of the parallelepiped, in crystallographic coordinates.
% {x1}: the opposite vertex.
% {nmol}: if nmol is positive, then detect the connected molecules in the 
%         selected region using mol_burst and keep only the molecules with
%         nmol or mol atoms.
%
% Required output variables:
% {mol}: the molecular description.
%

bohr2angstrom = 0.52917720859;

## crystal to cartesian
if (isfield(cr,"r"))
  r = cr.r;
else
  if (isfield(cr,"g"))
    g = cr.g;
  else
    cc = cos(cr.b);
    g = cr.a' * cr.a;
    g(1,2) = g(2,1) = g(1,2) * cc(3);
    g(1,3) = g(3,1) = g(1,3) * cc(2);
    g(2,3) = g(3,2) = g(2,3) * cc(1);
  endif
  r = chol(g)';
endif

if (LOG>0)
   printf("Crys2Cart matrix:\n");
   for i = 1 : 3
      for j = 1 : 3
         printf(" %15.9f", r(i,j));
      endfor
      printf("\n");
   endfor
endif

## find the relevant lattice vectors
ix0 = floor(x0); ix1 = ceil(x1)-1;
nvecs = prod(ix1-ix0+1); nvecsm1 = max(prod(ix1-1-(ix0+1)+1),0);
lvecs = zeros(nvecs,3);
n = 0;
for ix = ix0(1):ix1(1)
  for iy = ix0(2):ix1(2)
    for iz = ix0(3):ix1(3)
      n += 1;
      lvecs(n,:) = [ix iy iz];
    endfor
  endfor
endfor

## build the molecule, assume that cells on the border contribute 1/4
## of the atoms in the unit cell to save memory.
mol = struct();
mol.atname = cell();
mol.atnumber = zeros(1,nvecsm1*cr.nat+nvecs*cr.nat/4);
mol.atxyz = zeros(3,nvecsm1*cr.nat+nvecs*cr.nat/4);
n = 0;
for i = 1:cr.nat
  for j = 1:nvecs
    x = cr.x(i,:) + lvecs(j,:);
    if (all(x > x0) && all(x < x1))
      n += 1;
      mol.atname{n} = cr.attyp{cr.typ(i)};
      mol.atnumber(n) = cr.ztyp(cr.typ(i));
      mol.atxyz(:,n) = x * r;
    endif
  endfor
endfor
if (n > 0)
  mol.atnumber = mol.atnumber(1:n);
  mol.atxyz = mol.atxyz(:,1:n) * bohr2angstrom;
else
  mol.atnumber = mol.atxyz = [];
endif

if (nmol > 0)
  smol = mol_burst(mol);
  mol = molecule();
  for i = 1:length(smol)
    if (length(smol{i}.atnumber) >= nmol)
      mol = mol_merge(mol,smol{i});
    endif
  endfor
endif

if (isfield(cr,"name") && !isempty(cr.name))
  mol.name = cr.name;
endif

if (LOG>0)
  printf ("CRYSTALBOX (atoms: %d)\n",n);
  printf("-ind- -typ- --Z--            ---cart-coord---\n");
  for i = 1:n
    printf ("%5d %5s %5d",i, mol.atname{i},mol.atnumber(i));
    printf (" %12.6f %12.6f %12.6f\n", mol.atxyz(1:3,i));
  endfor
endif

endfunction
