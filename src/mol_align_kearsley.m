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

function [molout u rmsd] = mol_align_kearsley(mol0, mol1)
 % function [molout u rmsd] = mol_align_kearsley(mol0, mol1)
 %
 % mol_align_kearsley - rotate the mol0 molecule in order to minimize
 % the rmsd of the atomic positions compared to the mol1 molecule. The
 % two molecules must have the same number of atoms and the atoms must
 % be in the same order. The output molecule is the mol0 molecule,
 % rotated and translated to the mol1 center of mass.  The rotation
 % matrix u, and the RMSD of the atomic positions (rmsd) are also
% given.  This routine uses Kearsley's algorithm. (S.K. Kearsley, Acta
% Cryst. A 45 (1989) 208--210.)
%
% Input variables:
% mol0: rotated molecule
% mol1: target molecule.
% allow_inversion: if 1, permit the application of an inversion
% operation to match the two molecules. Warning: if mol0 is chiral, it
% is converted to its enantiomer.
%
% Output variables:
% molout: rotated and translated mol0.
% u: rotation matrix.
% rmsd: 
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: December 2011

  if (mol0.nat != mol1.nat)
    error("inconsistent number of atoms in mol0 and mol1")
  endif
  n = mol1.nat;
  xcm0 = mol_cmass(mol0);
  xcm1 = mol_cmass(mol1);
  mol0.atxyz = mol0.atxyz - xcm0 * ones(1,mol0.nat);
  mol1.atxyz = mol1.atxyz - xcm1 * ones(1,mol1.nat);

  ## xm and xp vectors:
  xm = mol0.atxyz - mol1.atxyz;
  xp = mol0.atxyz + mol1.atxyz;

  ## s matrix
  S = zeros(4,4);
  for i = 1:n
    S(1,1) = S(1,1) + (xm(1,i)^2 + xm(2,i)^2 + xm(3,i)^2);
    S(1,2) = S(1,2) + (xp(2,i)*xm(3,i) - xm(2,i)*xp(3,i));
    S(1,3) = S(1,3) + (xm(1,i)*xp(3,i) - xp(1,i)*xm(3,i));
    S(1,4) = S(1,4) + (xp(1,i)*xm(2,i) - xm(1,i)*xp(2,i));
    S(2,2) = S(2,2) + (xp(2,i)^2 + xp(3,i)^2 + xm(1,i)^2);
    S(2,3) = S(2,3) + (xm(1,i)*xm(2,i)-xp(1,i)*xp(2,i));
    S(2,4) = S(2,4) + (xm(1,i)*xm(3,i)-xp(1,i)*xp(3,i));
    S(3,3) = S(3,3) + (xp(1,i)^2+xp(3,i)^2+xm(2,i)^2);
    S(3,4) = S(3,4) + (xm(2,i)*xm(3,i)-xp(2,i)*xp(3,i));
    S(4,4) = S(4,4) + (xp(1,i)^2+xp(2,i)^2+xm(3,i)^2);
  endfor
  S(2,1) = S(1,2);
  S(3,1) = S(1,3);
  S(3,2) = S(2,3);
  S(4,1) = S(1,4);
  S(4,2) = S(2,4);
  S(4,3) = S(3,4);

  ## Diagonalize S and get the quaternions for the transformation:
  [vec,val] = eig(S);
  [sval,isort] = sort(diag(val));
  q = vec(1:4,isort(1));

  ## Convert the quaternion into the rotation matrix:
  u = zeros(3,3);
  u(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
  u(1,2) = 2*(q(2)*q(3) + q(1)*q(4));
  u(1,3) = 2*(q(2)*q(4) - q(1)*q(3));
  u(2,1) = 2*(q(2)*q(3) - q(1)*q(4));
  u(2,2) = q(1)^2 + q(3)^2 - q(2)^2 - q(4)^2;
  u(2,3) = 2*(q(3)*q(4) + q(1)*q(2));
  u(3,1) = 2*(q(2)*q(4) + q(1)*q(3));
  u(3,2) = 2*(q(3)*q(4) - q(1)*q(2));
  u(3,3) = q(1)^2 + q(4)^2 - q(2)^2 - q(3)^2;

  ## output quantities
  atxyz0 = u * mol0.atxyz;
  rmsd = sqrt(sum(sum((atxyz0 - mol1.atxyz).^2)) / mol0.nat);
  molout = mol0;
  molout.atxyz = atxyz0 + xcm1 * ones(1,n);

endfunction
