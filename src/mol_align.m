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

function [mol2new,R] = mol_align(mol1, iat1, mol2, iat2, LOG=1)
% function [mol2new,R] = mol_align(mol1, iat1, mol2, iat2, LOG=1)
%
% mol_align - translate and rotate mol2 for best aligment with mol1.
%
% The algorithm is based on:
% S.K. Kearsley, "On the orthogonal transformation used for structural
% comparisons", Acta Cryst. A 45 (1989) 208--210.
%
% Required input variables:
% {mol1,mol2}: the two molecules to align.
% {iat1,iat2}: list of the atom positions in the two molecules that should
%     align. This can be used to select part of the molecules or adapt
%     different orders of their atoms. It is required that iat1 and iat2
%     be of the same size. If iat1 and iat2 are empty, the whole molecules
%     will be used for alignment, and mol1 and mol2 must have equal number
%     of atoms.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {mol2new}: mol2 with the new alignment.
% {R}: rotation matrix to align mol2 with mol1. Notice that this *must*
%      be applied after the center of mass of mol2 (cm2) has been placed
%      at the origin. In other words, the sequence should be: translate
%      cm2 to the origin, then apply R, then move cm2 to its final
%      position.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

# Check input data, particularly iat1 and iat2.
if (isempty(iat1) && isempty(iat2))
   n1 = mol1.nat;
   n2 = mol2.nat;
   if (n1 == n2)
      iat1 = 1:n1; iat2 = iat1;
   else
      error('mol_align: mol sizes do nor match!');
   endif
elseif (!isequal(size(iat1), size(iat2)))
   error('mol_align: iat sizes do nor match!');
else
   n1 = length(iat1); n2 = length(iat2);
endif
if (n1 < 3)
   error('mol_align: at least three non-coplanar atoms are needed!');
endif

# Get coordinates that must be compared:
x1 = mol1.atxyz(1:3,iat1); mass1 = mol1.atmass(iat1);
x2 = mol2.atxyz(1:3,iat2); mass2 = mol2.atmass(iat2);

# Get the centers of mass and put them at the origin:
cm1 = x1 * mass1' / sum(mass1);
cm2 = x2 * mass2' / sum(mass2);
x1 = x1 - cm1 * ones(1,n1);
x2 = x2 - cm2 * ones(1,n2);


## # Average the non-singular matrices that convert triplets of mol2 atoms
## # into their equivalent in mol1:
## R = zeros(3,3);
## m = 0;
## for i = 1 : n1
##    for j = i+1 : n1
##       for k = j+1 : n1
##          ijk = [i,j,k];
##          B = x2(1:3,ijk);
##          if (abs(det(B)) > 1e-3)
##             m++;
##             R = R + x1(1:3,ijk) * inv(B);
##          endif
##       endfor
##    endfor
## endfor
## if (m > 0)
##    R = R / m;
## else
##    error('mol_align: singular transform matrix!');
## endif

# xm and xp vectors:
xm = x1-x2;
xp = x1+x2;
S = zeros(4,4);
for i = 1 : n1
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

# Diagonalize S and get the quaternions for the transformation:
[vec,val] = eig(S);
[sval,isort] = sort(diag(val));
q = vec(1:4,isort(1));

# Convert the quaternion into the rotation matrix:

R = zeros(3,3);
R(1,1) = q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2;
R(1,2) = 2*(q(2)*q(3) + q(1)*q(4));
R(1,3) = 2*(q(2)*q(4) - q(1)*q(3));
R(2,1) = 2*(q(2)*q(3) - q(1)*q(4));
R(2,2) = q(1)^2 + q(3)^2 - q(2)^2 - q(4)^2;
R(2,3) = 2*(q(3)*q(4) + q(1)*q(2));
R(3,1) = 2*(q(2)*q(4) + q(1)*q(3));
R(3,2) = 2*(q(3)*q(4) - q(1)*q(2));
R(3,3) = q(1)^2 + q(4)^2 - q(2)^2 - q(3)^2;

# Transform the whole mol2:
mol2new = mol2;
m2 = mol2.nat;
mol2new.atxyz = R * (mol2.atxyz - cm2 * ones(1,m2)) + cm1 * ones(1,m2);

if (LOG > 0)
   printf('mol_align: alignment of fragments of two molecules\n');
   printf('Fragment %d (%d):', 1, mol1.nat);
   printf(' %d', iat1);
   printf('\n');
   printf('Fragment %d (%d):', 2, mol2.nat);
   printf(' %d', iat2);
   printf('\n');
   printf('Quaternion: %.9f %.9f %.9f %.9f\n', q);
   printf('Rotation matrix:\n');
   printf('   %15.9f %15.9f %15.9f\n', R);
   printf('\nPattern molecule:\n');
   for i = 1 : n1
      ii = iat1(i);
      printf('%4d %2s', ii, mol1.atname{ii});
      printf('%12.6f%12.6f%12.6f', mol1.atxyz(1:3,ii));
      printf('\n');
   endfor
   printf('\nMolecule to align:\n');
   for i = 1 : n1
      ii = iat2(i);
      printf('%4d %2s', ii, mol2.atname{ii});
      printf('%12.6f%12.6f%12.6f', mol2.atxyz(1:3,ii));
      printf('%12.6f%12.6f%12.6f', mol2new.atxyz(1:3,ii));
      printf('\n');
   endfor
endif

endfunction
