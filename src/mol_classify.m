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

function [iclass] = mol_classify(mol, iat, bondfactor=1.2, LOG=1)
% function [iclass] = mol_classify(mol, iat, bondfactor=1.2, LOG=1)
%
% mol_classify - Group into classes the atoms of mol in the iat list.
%
% Classification is based on:
%  1. atomic number. Ordered heavier to lighter.
%  2. equivalence of binding environments.
%
% Required input variables:
% {mol}: molecule to classify.
% {iat}: list of the atoms included inthe classification.
%     If iat is empty, the whole molecule will be used.
%
% Optional input variables (all have default values):
% {bondfactor=1.20}: Extra allowance to consider two atoms bonded.
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {iclass}: class for each atom. Atoms with iclass<0 have not been classed.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2012

   global dbdefined
   global atdb

# Check iat data.
nat = length(mol.atnumber);
iclass = -ones(1,nat);
if (isempty(iat))
   is = 1:nat;
else
   is = iat;
endif
nis = length(is);
for i = 1:nis
   iclass(is(i)) = 0;
endfor

# Classify atoms according to
#  1. atomic number. Ordered heavier to lighter.
#  2. equivalence of binding environments.

# Sort by atomic order:
for i = 1:nis
   for j = i+1:nis
      if (mol.atnumber(is(j)) > mol.atnumber(is(i)))
         ii = is(i);
         is(i) = is(j);
         is(j) = ii;
      endif
   endfor
endfor

# Get connection matrix:
# (careful: we keep the natural order for the atoms).
conn = zeros(nis,nis);
for i = 1:nis
   dist(i,i) = 0.0;
   zi = z(i) = mol.atnumber(i);
   for j = 1:i-1
      x = mol.atxyz(1:3,i) - mol.atxyz(1:3,j);
      dist(i,j) = dist(j,i) = sqrt(x' * x);
      zj = mol.atnumber(j);
      dconn = bondfactor * (atdb.rcov(zi) + atdb.rcov(zj));
      conn(i,j) = conn(j,i) = (dist(i,j) <= dconn);
   endfor
endfor

# Get bond index:
# bn() number of bonds associated to an atom
# bi() the weight of a bond is the tomic number of the bonded atom
for i = 1:nis
   i1 = is(i);
   bn(i) = ones(1,nis) * conn(:,i);
   bi(i) = z(:)' * conn(:,i);
   %printf("%2d (%d) %d - %d\n", i, z(i), bn(i), bi(i));
   %for j = 1:nis; printf("%3d",z(j)); endfor; printf("\n");
   %for j = 1:nis; printf("%3d",conn(i,j)); endfor; printf("\n");
   %printf("\n");
endfor

# Sort by bi's as a second criterion:
istart = 1;
zstart = z(istart);
iend = istart;
for i = 2:nis
   i1 = is(i);
   zi = z(i1);
   if (zi == zstart && i < nis)
      iend = i;
   else
      %printf("[%d,%d] sort by bi\n", istart, iend);
      for k1 = istart:iend
         kk1 = is(k1);
         for k2 = istart:k1-1
            kk2 = is(k2);
            if (bi(kk2) < bi(kk1))
               ii = is(k1);
               is(k1) = is(k2);
               is(k2) = ii;
            endif
         endfor
      endfor
      istart = iend = i;
      zstart = zi;
   endif
endfor

if (LOG > 0)
   printf('mol_classify: sort atoms in classes\n');
   for i = 1:nat
      i1 = is(i);
      z1 = mol.atnumber(i1);
      s1 = mol.atname{i1};
      b1 = bi(i1);
      printf("%4d - %4d %3d %3s %3d\n", i, i1, z1, s1, b1);
   endfor
   printf('\n');
endif

endfunction
