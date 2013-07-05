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

function [mol] = cr_molmotif(cr, allinmaincell=1, ilist=[], bondfactor=1.20, LOG=0)
% function [mol] = cr_molmotif(cr, allinmaincell=1, ilist=[], bondfactor=1.20, LOG=1)
%
% cr_molmotif - return one or several molecular motifs where, at least,
% one atom is within the main unit cell.
%
% Required input variables:
% {cr}: struct containing the crystal.
%
% Optional input variables (all have default values):
% {allinmaincell}: determine *all* the molecular motifs that have, at least,
%     an atom in the main cell. Default if ilist is empty.
% {ilist[]}: list of atoms that will be used to determine molecular motifs.
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {mol}: description of the molecular motifs, including the cartesian
%      coordinates of the atoms.
%

bohr2angstrom = 0.52917720859;

crys2car = zeros(3,3);
ared = cr.a(1); bred = cr.a(2); cred = cr.a(3);
calph = cos(cr.b(1)); cbeta = cos(cr.b(2)); cgamm = cos(cr.b(3));
HH1 = sqrt(1d0-cgamm*cgamm);
HH2 = (calph-cbeta*cgamm) / sqrt(1d0-cgamm*cgamm);
HH3 = sqrt(1d0 - cbeta*cbeta - HH2*HH2);
crys2car(1,1) = ared;
crys2car(1,2) = bred * cgamm;
crys2car(1,3) = cred * cbeta;
crys2car(2,2) = bred * HH1;
crys2car(2,3) = cred * HH2;
crys2car(3,3) = cred * HH3;
crys2car = crys2car * bohr2angstrom;
if (LOG>0)
   printf("Crys2Cart matrix:\n");
   for i = 1 : 3
      for j = 1 : 3
         printf(" %15.9f", crys2car(i,j));
      endfor
      printf("\n");
   endfor
endif

# Cell lengths are in bohr, cartesians should be in angstrom

# Atoms in the main cell
# natoms = rows(cr.x);
natoms = cr.nat;
ZZat = cr.ztyp(cr.typ);
Symbat = cr.attyp(cr.typ);

# Get atoms in the main cell and the neigbohr cells:
nall = 0;
step = 1;
if (allinmaincell == 1)
   inmotif = zeros(1,27*natoms);
   istep = zeros(1,27*natoms);
   for i = [0, 1, -1]
      for j = [0, 1, -1]
         for k = [0, 1, -1]
            main = (i==0 && j==0 && k==0);
            for l = 1 : natoms
               nall++;
               if (main)
                  inmotif(nall) = 1;
                  istep(nall) = step;
               endif
               atindex(nall) = l;
               xyz(1:3,nall) = cr.x(l,1:3)' + [i;j;k];
            endfor
         endfor
      endfor
   endfor
elseif (isempty(ilist))
   inmotif = zeros(1,27*natoms);
   istep = zeros(1,27*natoms);
   for i = [0, 1, -1]
      for j = [0, 1, -1]
         for k = [0, 1, -1]
            main = (i==0 && j==0 && k==0);
            for l = 1 : natoms
               nall++;
               if (main && l == 1)
                  inmotif(nall) = 1;
                  istep(nall) = step;
               endif
               atindex(nall) = l;
               xyz(1:3,nall) = cr.x(l,1:3)' + [i;j;k];
            endfor
         endfor
      endfor
   endfor
else
   inmotif = zeros(1,26*natoms+length(ilist));
   istep = zeros(1,27*natoms);
   for i = [0, 1, -1]
      for j = [0, 1, -1]
         for k = [0, 1, -1]
            main = (i==0 && j==0 && k==0);
            if (main)
               for l = ilist
                  nall++;
                  inmotif(nall) = 1;
                  istep(nall) = step;
                  atindex(nall) = l;
                  xyz(1:3,nall) = cr.x(l,1:3)' + [i;j;k];
               endfor
            else
               for l = 1 : natoms
                  nall++;
                  inmotif(nall) = 0;
                  istep(nall) = 0;
                  atindex(nall) = l;
                  xyz(1:3,nall) = cr.x(l,1:3)' + [i;j;k];
               endfor
            endif
         endfor
      endfor
   endfor
endif
xyz(1:3,1:nall) = crys2car(1:3,1:3) * xyz(1:3,1:nall);
printf("DBG: nall %d [%d,%d]\n", nall, size(xyz));

for i = 1:nall
   i1 = atindex(i);
   zi = ZZat(i1);
   [symb,atom] = mol_dbsymbol(zi);
   radius(i) = atom.rcov;
   %%% printf("DBG at(%5d,%2s): %7.2f\n", i, symb, radius(i));
endfor

# Do iteratively until all seed atoms have been used
do
   new = 0;
   step++;
   listin = find(inmotif>=1);
   listout = find(inmotif<=0);
   printf("DBG: step(%d) - in: %d %d\n", step, length(listin), length(listout));
   for i = listin
      for j = listout
         if (i > nall || j > nall)
            printf("DBG: %d %d [%d,%d]\n", i, j, size(xyz));
         endif
         xx = xyz(1:3,i) - xyz(1:3,j);
         dij = sqrt(xx' * xx);
         bonded = dij <= bondfactor*(radius(i)+radius(j));
       % if (bonded)
       %    ss = "Yes";
       % else
       %    ss = "N";
       % endif
       % printf("DBG: d(%4d,%4d): %12.4f (%6.2f,%6.2f) %s\n",i,j,dij,radius(i),radius(j),ss);
         if (bonded)
            new++;
            inmotif(j) = 1;
            istep(j) = step;
         endif
      endfor
   endfor
until (new <=0)

if (LOG>0)
   listin = find(inmotif>=1);
   printf ("MOLMOTIF (bond factor: %.3f, atoms: %d)\n", bondfactor, length(listin));
   printf ("-ord- -ind- --Z-- -stp- ---cart-coord--...\n");
   n = 0;
   for i = listin
      i1 = atindex(i);
      zi = ZZat(i1);
      printf ("%5d %5d %5d %5d", ++n, i, zi, istep(i));
      printf (" %15.5f %15.5f %15.5f\n", xyz(1:3,i));
   endfor
   cm = sum(xyz(1:3,listin),2)/n;
   printf ("Center of Mass: %15.5f %15.5f %15.5f\n", cm);
endif

# Create the final mol output:
listin = find(inmotif>=1);
n = 0;
for i = listin
   ++n;
   i1 = atindex(i);
   zi = ZZat(i1);
   [symb,atom] = mol_dbsymbol(zi);
   mol.atname{n} = symb;
   mol.atnumber(n) = zi;
   mol.atmass(n) = atom.mass;
   mol.atxyz(1:3,n) = xyz(1:3,i);
endfor

endfunction
