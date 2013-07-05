% Copyright (C) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [mol2new,R] = mol_opt_align(mol1, iat1, mol2, iat2, bondfactor=1.2, LOG=1)
% function [mol2new,R] = mol_opt_align(mol1, iat1, mol2, iat2, bondfactor=1.2, LOG=1)
%
% mol_opt_align - translate and rotate mol2 for best aligment with mol1.
%
% The align algorithm is based on:
% S.K. Kearsley, "On the orthogonal transformation used for structural
% comparisons", Acta Cryst. A 45 (1989) 208--210.
% This requires that equivalent atoms within the two molecules be
% identified by hand, an easy task if done for a reduced number of
% molecules with a limited number of atoms. However, when done for a larger
% number of molecules and atoms, an automated algorithm is most welcome. This
% is the purpose of the present routine.
%
% Required input variables:
% {mol1,mol2}: the two molecules to align.
% {iat1,iat2}: list of the atom positions in the two molecules that should
%     align. This is only used to select parts of the molecules. To align
%     molecules with the strict orders of iat1 and iat2 use the mol_align()
%     routine.
%     It is required that iat1 and iat2 be of the same size. If iat1 and
%     iat2 are empty, the whole molecules will be used for alignment,
%     and mol1 and mol2 must have equal number of atoms.
%
% Optional input variables (all have default values):
% {bondfactor=1.20}: Extra allowance to consider two atoms bonded.
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
% Created: December 2012

   global dbdefined
   global atdb

# Check input data, particularly iat1 and iat2.
if (isempty(iat1) && isempty(iat2))
   n1 = length(mol1.atnumber);
   n2 = length(mol2.atnumber);
   if (n1 == n2)
      iat1 = 1:n1; iat2 = iat1;
   else
      error('mol_opt_align: mol sizes do nor match!');
   endif
elseif (!isequal(size(iat1), size(iat2)))
   error('mol_opt_align: iat sizes do nor match!');
else
   n1 = length(iat1); n2 = length(iat2);
endif
if (n1 < 3)
   error('mol_opt_align: at least three non-linear atoms are needed!');
endif

# Get coordinates that must be compared:
x1 = mol1.atxyz(1:3,iat1); mass1 = mol1.atmass(iat1);
x2 = mol2.atxyz(1:3,iat2); mass2 = mol2.atmass(iat2);

# Get the centers of mass and put them at the origin:
cm1 = x1 * mass1' / sum(mass1);
cm2 = x2 * mass2' / sum(mass2);
x1 = x1 - cm1 * ones(1,n1);
x2 = x2 - cm2 * ones(1,n2);

# Classify atoms in both lists.
# The sequence of criteria is:
#  1. atomic number. Ordered heavier to lighter.
#  2. equivalence of binding environments.

is1 = iat1;
for i = 1:n1
   for j = i+1:n1
      if (mol1.atnumber(is1(j)) > mol1.atnumber(is1(i)))
         ii = is1(i);
         is1(i) = is1(j);
         is1(j) = ii;
      endif
   endfor
endfor

b1 = zeros(n1,n1);
for i = 1:n1
   dist(i,i) = 0.0;
   zi = mol1.atnumber(i);
   for j = 1:i-1
      x = mol1.atxyz(1:3,i) - mol1.atxyz(1:3,j);
      dist(i,j) = dist(j,i) = sqrt(x' * x);
      zj = mol1.atnumber(j);
      dconn = bondfactor * (atdb.rcov(zi) + atdb.rcov(zj));
      conn(i,j) = conn(j,i) = (dist(i,j) <= dconn);
   endfor
endfor
mol1.conn = conn;

class = 1; iclass(class) = is1(1); zc = mol1.atnumber(is1(i));
ic1(1) = class;
for i = 2 : n1
   if (mol1.atnumber(is1(i)) == zc)
      ic1(i) = class;
   else
      iclass(++class) = is1(i);
      zc = mol1.atnumber(is1(i));
      ic1(i) = class;
   endif
endfor

is2 = iat2;
for i = 1:n2
   for j = i+1:n2
      if (mol2.atnumber(is2(j)) > mol2.atnumber(is2(i)))
         ii = is2(i);
         is2(i) = is2(j);
         is2(j) = ii;
      endif
   endfor
endfor

b1 = zeros(n2,n2);
for i = 1:n2
   dist(i,i) = 0.0;
   zi = mol2.atnumber(i);
   for j = 1:i-1
      x = mol2.atxyz(1:3,i) - mol2.atxyz(1:3,j);
      dist(i,j) = dist(j,i) = sqrt(x' * x);
      zj = mol2.atnumber(j);
      dconn = bondfactor * (atdb.rcov(zi) + atdb.rcov(zj));
      conn(i,j) = conn(j,i) = (dist(i,j) <= dconn);
   endfor
endfor
mol2.conn = conn;

class = 1; iclass(class) = is2(1); zc = mol2.atnumber(is2(i));
ic2(1) = class;
for i = 2 : n1
   if (mol2.atnumber(is2(i)) == zc)
      ic2(i) = class;
   else
      iclass(++class) = is2(i);
      zc = mol1.atnumber(is2(i));
      ic2(i) = class;
   endif
endfor

is2 = iat2;
for i = 1:n2
   for j = i+1:n2
      if (mol2.atnumber(is2(j)) > mol2.atnumber(is2(i)))
         ii = is2(i);
         is2(i) = is2(j);
         is2(j) = ii;
      endif
   endfor
endfor

if (LOG > 0)
   printf('mol_align_opt: alignment of fragments of two molecules\n');
   for i = 1:n1
      i1 = is1(i); i2 = is2(i);
      ###i1 = ic1(i); i2 = ic2(i);
      z1 = mol1.atnumber(i1); z2 = mol2.atnumber(i2);
      s1 = mol1.atname{i1}; s2 = mol2.atname{i2};
      printf("%4d - %4d %3d %3s - %4d %3d %3s\n", i, i1, z1, s1, i2, z2, s2);
   endfor
   printf('\n');
endif

endfunction
