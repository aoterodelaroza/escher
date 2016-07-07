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

function [mol] = mol_internalgeometry (mol, bondfactor=1.15, LOG=2)
% function [mol] = mol_internalgeometry (mol,bondfactor=1.15,LOG=1)
%
% mol_internalgeometry - determines the internal geometry of the molecule
% or molecular fragmente, i.e. the atoms bonded, bond distances, angles
% and dihedral angles. The algorithm assumes two atoms to be bonded if
% their distance is smaller or equal to the sum of the atomic radii
% multiplied by a correction factor (the "bondfactor" parameter).
%
% Required input variables:
% mol: molecular database including the cartesian coordinates of all atoms.
%
% Optional input variables (all have default values):
% {LOG = 1}: print results if LOG>0.
%    LOG = 0  no output.
%    LOG = 1  distances and connectivity matrix.
%    LOG = 2  angles, linear/planar molecules and rings
%    LOG = 3  dihedral angles (very expensive if there are many atoms)
% {bondfactor=1.15}: Extra allowance to consider two atoms bonded.
%
% Output:
% mol: molecular database with the internal geometry added.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: July 2011

   global dbdefined
   global atdb
   bohrtoans = 0.52917720859;
#   DEBUG = 1;

   # Compute the matrix of distances and the connectivity matrix:
   M = mol.nat;
   dist = conn = zeros(M);
   for i = 1:M
      dist(i,:) = norm(mol.atxyz - mol.atxyz(:,i) * ones(1,M),2,"columns");
   endfor
   for i = 1:M-1
      zi = mol.atnumber(i);
      conn(i,i) = 0;
      for j = i+1:M
         zj = mol.atnumber(j);
         dconn = bondfactor * (atdb.rcov(zi) + atdb.rcov(zj));
         conn(i,j) = conn(j,i) = (dist(i,j)<=dconn);
#         if (DEBUG>0)
#            printf('%s-%d,%s-%d ', mol.atname{i}, i, mol.atname{j}, j);
#            printf('%10.5f %10.5f ', dist(i,j), dconn);
#            printf('conn:%d\n', conn(i,j));
#         endif
      endfor
   endfor
   mol.conn = conn;

   # Check for unconnected fragments in the molecule:
   taken = zeros(1,M);
   nfrag = 0;
   for i = 1:M
      if (taken(i) == 0)
         taken(i) = 1;
         nfrag = nfrag + 1;
         frag.nat(nfrag) = 1;
         frag.iat{nfrag}(1) = i;
         do
            inew = 0;
            for j = i+1 : M
               if (taken(j) == 0)
                  isinfrag = 0;
                  for k = 1 : frag.nat(nfrag)
                     k1 = frag.iat{nfrag}(k);
                     if (conn(j,k1) > 0)
                        isinfrag = 1;
                     endif
                  endfor
                  if (isinfrag > 0)
                     nf = ++frag.nat(nfrag);
                     frag.iat{nfrag}(nf) = j;
                     taken(j) = 1;
                     inew = inew + 1;
                  endif
               endif
            endfor
         until (inew == 0)
      endif
   endfor
   mol.nfrag = nfrag;
   mol.ifrag = frag.iat;

   # Determine cm of each fragment:
   for i = 1 : nfrag
      x = mol.atxyz(1:3,frag.iat{i});
      mass = mol.atmass(frag.iat{i});
      cm = x * mass' / sum(mass);
      frag.cm{i} = cm;
   endfor
   mol.cmfrag = frag.cm;

   if (LOG > 0)
      printf("(Coordinates in Angstrom)\n");
      printf("    Atom  number  cn       x              y              z\n");
      for i = 1:M
         ni = cell2mat(mol.atname(i));
         cni = sum(conn(i,:));
         zi = mol.atnumber(i);
         x = mol.atxyz(1:3,i);
         printf("%3d %-6s%4d%4d  %15.10f%15.10f%15.10f\n", i, ni, zi, cni, x);
      endfor

      #
      # Print distances between bonded atoms:
      printf("\nBond distances (Angstrom):\n");
      for i = 1:M
         si = sprintf("%s(%d)", cell2mat(mol.atname(i)), i);
         for j = 1:i-1
            if (conn(i,j))
               sj = sprintf("%s(%d)", cell2mat(mol.atname(j)), j);
               printf("BOND %12.6f %8s-%s\n", dist(i,j), si, sj);
            endif
         endfor
      endfor
   endif

   if (LOG > 1)
      #
      # Print angles between bonds:
      printf("\nBond angles (degrees):\n");
      for i = 1:M
         for j = 1:i-1
            for k = 1:j-1
               if (conn(i,j) & conn(j,k))
                  t = [i, j, k];
               elseif (conn(i,k) & conn(k,j))
                  t = [i, k, j];
               elseif (conn(j,i) & conn(i,k))
                  t = [j, i, k];
               else
                  continue
               endif
               x21 = mol.atxyz(1:3,t(2)) - mol.atxyz(1:3,t(1));
               x23 = mol.atxyz(1:3,t(2)) - mol.atxyz(1:3,t(3));
               pesc = dot(x21,x23);
               pvec = cross(x21,x23);
               pvec = sqrt(dot(pvec,pvec));
               angle = atan2(pvec,pesc)*180/pi;
               s1 = sprintf("%s(%d)",cell2mat(mol.atname(t(1))),t(1));
               s2 = sprintf("%s(%d)",cell2mat(mol.atname(t(2))),t(2));
               s3 = sprintf("%s(%d)",cell2mat(mol.atname(t(3))),t(3));
               printf("ANGLE %12.6f %s-%s-%s\n", angle, s1, s2, s3);
            endfor
         endfor
      endfor
   endif

   if (LOG > 3)
      #
      # Print dihedral angles between bonds:
      printf("\nDihedral angles (degrees):\n");
      for i = 1:M
         for j = 1:i-1
            for k = 1:j-1
               for l = 1:k-1
                  if (conn(i,j) & conn(j,k) & conn(k,l))
                     t = [i, j, k, l];
                  elseif (conn(j,i) & conn(i,k) & conn(k,l))
                     t = [j, i, k, l];
                  elseif (conn(i,k) & conn(k,j) & conn(j,l))
                     t = [i, k, j, l];
                  elseif (conn(j,k) & conn(k,i) & conn(i,l))
                     t = [j, k, i, l];
                  elseif (conn(k,i) & conn(i,j) & conn(j,l))
                     t = [k, i, j, l];
                  elseif (conn(k,j) & conn(j,i) & conn(i,l))
                     t = [k, j, i, l];
                  elseif (conn(i,j) & conn(j,l) & conn(l,k))
                     t = [i, j, l, k];
                  elseif (conn(j,i) & conn(i,l) & conn(l,k))
                     t = [j, i, l, k];
                  elseif (conn(i,k) & conn(k,l) & conn(l,j))
                     t = [i, k, l, j];
                  elseif (conn(j,k) & conn(k,l) & conn(l,i))
                     t = [j, k, l, i];
                  elseif (conn(k,i) & conn(i,l) & conn(l,j))
                     t = [k, i, l, j];
                  elseif (conn(k,j) & conn(j,l) & conn(l,i))
                     t = [k, j, l, i];
                  else
                     continue
                  endif
                  x12 = mol.atxyz(1:3,t(2)) - mol.atxyz(1:3,t(1));
                  x23 = mol.atxyz(1:3,t(3)) - mol.atxyz(1:3,t(2));
                  x34 = mol.atxyz(1:3,t(4)) - mol.atxyz(1:3,t(3));
                  n123 = cross(x12,x23);
                  n234 = cross(x23,x34);
                  pmod = sqrt(dot(n123,n123) * dot(n234,n234));
                  if (pmod > 0)
                     pesc = dot(n123,n234);
                     pvec = cross(n123,n234);
                     pvec = sqrt(dot(pvec,pvec));
                     dihedr = atan2(pvec,pesc) * 180/pi;
                  else
                     dihedr = -720 ### Error!!!
                  endif
                  s1 = sprintf("%s(%d)",cell2mat(mol.atname(t(1))),t(1));
                  s2 = sprintf("%s-%s(%d)",s1,cell2mat(mol.atname(t(2))),t(2));
                  s3 = sprintf("%s-%s(%d)",s2,cell2mat(mol.atname(t(3))),t(3));
                  s4 = sprintf("%s-%s(%d)",s3,cell2mat(mol.atname(t(4))),t(4));
                  printf("DIHEDRAL %12.6f %s\n", dihedr, s4);
               endfor
            endfor
         endfor
      endfor
   endif

   if (LOG > 0)
      #
      # Connection matrix:

      #
      # Number of unconnected fragments:
      printf("\n");
      for i = 1 : nfrag
         printf("Fragment (%02d):", i);
         printf(" %d", frag.iat{i});
         printf("\n");
         printf("frag. center of mass: ");
         printf("%15.10f %15.10f %15.10f", frag.cm{i});
         printf("\n");
      endfor
   endif

   if (LOG > 1)
      # Check if the molecule is linear or planar (or approximately):
      # --- xxx

      if (M > 100)
         printf("Identification of ring structures skipped due to the\n");
         printf("big number of atoms. The search could take ages!\n");
         return
      endif
      #
      # Identify rings of connected atoms:
      # Algorithmic limitation: Instead of searching for rings of arbitrary
      # size we will only examine the occurrence of n-rings, where n: 3--7.
      # The present algorithm can be improved and generalized
      # using recursivity.
      nring = 0;
      nsupring = 0;
      for i1 = 1:M
         for i2 = i1+1:M
            if (!conn(i1,i2)) continue endif
            for i3 = i2+1:M
               if (!conn(i2,i3)) continue endif
               if (conn(i2,i3) && conn(i3,i1))
                  # A 3-ring
                  ring{++nring} = [i1, i2, i3];
                  continue
               endif
               for i4 = i3+1:M
                  if (!conn(i3,i4)) continue endif
                  if (conn(i3,i4) && conn(i4,i1))
                     # A 4-ring
                     ring{++nring} = [i1, i2, i3, i4];
                     continue
                  endif
                  for i5 = i4+1:M
                     if (!conn(i4,i5)) continue endif
                     if (conn(i4,i5) && conn(i5,i1))
                        # A 5-ring
                        ring{++nring} = [i1, i2, i3, i4, i5];
                        continue
                     endif
                     for i6 = i5+1:M
                        if (!conn(i5,i6)) continue endif
                        if (conn(i5,i6) && conn(i6,i1))
                           # A 6-ring
                           ring{++nring} = [i1, i2, i3, i4, i5, i6];
                           continue
                        endif
                        for i7 = i6+1:M
                           if (!conn(i6,i7)) continue endif
                           if (conn(i6,i7) && conn(i7,i1))
                              # A 7-ring
                              ring{++nring} = [i1, i2, i3, i4, i5, i6, i7];
                              continue
                           endif
                           for i8 = i7+1:M
                              if (!conn(i7,i8)) continue endif
                              if (conn(i7,i8) && conn(i8,i1))
                                 # A 8-ring
                                 ring{++nring} = ...
                                     [i1, i2, i3, i4, i5, i6, i7, i8];
                                 continue
                              endif
                              # A possible ring with >8 atoms
                              # Unsupported here
                              nsupring++
                           endfor
                        endfor
                     endfor
                  endfor
               endfor
            endfor
         endfor
      endfor
      printf("\nPossible rings not determined: %6d\n", nsupring);
      printf("Rings found in the structure.: %6d\n", nring);
      for i = 1:nring
         nr = length(ring{i});
         printf("  %3d-> %d-ring formed by: ", i, nr);
         printf("%d ", ring{i});
         printf("\n");
         # Decompose the ring in succesive triangles and average the
         # normal vector from all the triangles:
         vn = [0;0;0];
         for k1 = 1:nr
            k2 = mod(k1,nr)+1;
            k3 = mod(k2,nr)+1;
            x12 = mol.atxyz(1:3,ring{i}(k2)) - mol.atxyz(1:3,ring{i}(k1));
            x23 = mol.atxyz(1:3,ring{i}(k3)) - mol.atxyz(1:3,ring{i}(k2));
            vn = vn + cross(x12,x23);
         endfor
         vmod = sqrt(dot(vn,vn));
         if (vmod > 0)
            vn = vn / vmod;
            printf("Average normal vector of the ring: ");
            printf("%12.6f %12.6f %12.6f\n", vn);
         else
            printf("Ring normal vector can't be normalized!\n");
         endif
      endfor
   endif


endfunction
