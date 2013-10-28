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

function mol = cr_spherebox(cr, x0, rad0, LOG=0)
% function mol = cr_spherebox(cr, x0, rad0, LOG=0)
%
% cr_spherebox - create a molecule (mol) by cutting the piece of the crystal
% contained within a sphere of radius r centered around x0 (cryst. coords.)
%
% Required input variables:
% {cr}: struct of the crystal.
% {x0}: center of the sphere (crystallographic coordinates).
% {rad0}: radius of the sphere (bohr). If an array [r0 r1] is provided,
%         restrict the selection to the region between the r0-sphere and
%         the r1-sphere.
%
% Optional input variables (all have default values):
% {LOG}: verbose flag.
%
% Required output variables:
% {mol}: the molecular description.
%

bohr2angstrom = 0.52917720859;
if (isscalar(rad0))
  rad = rad0;
  rmin = 0;
else
  rad = rad0(2);
  rmin = rad0(1);
endif

## transform x0 to the main cell
x0 = x0 - floor(x0);

## crystal to cartesian
r = cr.r;
g = cr.g;
rinv = inv(r);

## find the relevant lattice vectors
xc = x0 * r;
xc0 = xc - [rad rad rad];
xc1 = xc + [rad rad rad];
ix0 = abs(floor(xc0 * rinv));
ix1 = abs(ceil(xc1 * rinv));
ix0 = max(ix0,ix1)+1;

lvecs = zeros(prod(2*ix0+1),3);
n = 0;
for ix = -ix0(1):ix0(1)
  for iy = -ix0(2):ix0(2)
    for iz = -ix0(3):ix0(3)
      n += 1;
      lvecs(n,:) = [ix iy iz];
    endfor
  endfor
endfor
nvecs = size(lvecs,1);

## build the molecule, assume that cells on the border contribute 1/4
## of the atoms in the unit cell to save memory.
mol = molecule();
mol.atname = cell();
mol.atnumber = zeros(1,nvecs*cr.nat+nvecs*cr.nat/4);
mol.atxyz = zeros(3,nvecs*cr.nat+nvecs*cr.nat/4);
n = 0;
rmin2 = rmin * rmin; rad2 = rad * rad;
for i = 1:cr.nat
  for j = 1:nvecs
    x = cr.x(i,:) + lvecs(j,:);
    xd = x - x0;
    dist2 = xd * g * xd';
    if (dist2 <= rad2 && dist2 >= rmin2)
      n += 1;
      mol.atname{n} = cr.attyp{cr.typ(i)};
      mol.atnumber(n) = cr.ztyp(cr.typ(i));
      mol.atxyz(:,n) = x * r;
    endif
  endfor
endfor
mol.nat = n
if (n > 0)
  mol.atnumber = mol.atnumber(1:n);
  mol.atxyz = mol.atxyz(:,1:n) * bohr2angstrom;
else
  mol.atnumber = mol.atxyz = [];
endif
if (isfield(cr,"name") && !isempty(cr.name))
  mol.name = cr.name;
endif
mol = mol_fillatmass(mol);

endfunction
