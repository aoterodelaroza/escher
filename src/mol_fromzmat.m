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

function mol = mol_fromzmat(mol0, izmat, xzmat)
% function mol = mol_fromzmat(mol0, izmat, xzmat)
%
% mol_fromzmat - rebuild the atomic coordinates of the molecule from a
%   z-matrix.
%
% Required input variables:
% mol0: input molecule
% izmat: array of z-matrix indices
% xzmat: array of distances, angles, and dihedrals.
%
% Output:
% mol: output molecule.
%

  mol = mol0;

  ## prepare the array for the new atomic coordinates
  x = zeros(size(mol.atxyz))';

  ## first atom at zero
  if (mol.nat <= 1)
    return
  endif

  ## second atom along the z axis
  x(2,3) = xzmat(2,1);
  
  if (mol.nat <= 2)
    return
  endif

  ## third atom forming an angle with the other two, and in the xz plane
  x(3,:) = zmat_step(x(izmat(3,2),:),x(izmat(3,3),:),[1 0 0],xzmat(3,1),xzmat(3,2),0);

  ## rest of the atoms
  for i = 4:mol.nat
    x(i,:) = zmat_step(x(izmat(i,2),:),x(izmat(i,3),:),x(izmat(i,4),:),xzmat(i,1),xzmat(i,2),xzmat(i,3));
  endfor

  ## assign the new positions to the atoms in mol
  for i = 1:mol.nat
    mol.atxyz(:,izmat(i,1)) = x(i,:);
  endfor

endfunction
