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

function [] = mol_geom (mol, tlist, LOG=1)
% function [] = mol_geom (mol,tlist,LOG=1)
%
% mol_geom - given a molecule (mol) and a list of atoms this routine
% prints the distances, angles and dihedrals of the atoms connected as
% indicated in the list:
%
% Required input variables:
% mol: molecular database including the cartesian coordinates of all atoms.
% tlist: sequence of atoms in the molecule that form the desired connection.
%
% Optional input variables (all have default values):
% None
%
% Output:
% printed values of distances, angles, and dihedrals.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2012

   # Distances between pairs:
   nat = mol.nat;
   tat = length(tlist);
   if (tat > 1)
      printf("Distances between atoms:\n");
      for i = 1:tat-1
         t1 = tlist(i); t2 = tlist(i+1);
         x = mol.atxyz(1:3,t1) - mol.atxyz(1:3,t2);
         dd = sqrt(x' * x);
         s1 = sprintf("%s(%d)", cell2mat(mol.atname(t1)), t1);
         s2 = sprintf("%s(%d)", cell2mat(mol.atname(t2)), t2);
         printf("%12.6f %s-%s\n", dd, s1, s2);
      endfor
      printf("\n");
   endif
   if (tat > 2)
      printf("Angles between triplets:\n");
      for i = 1:tat-2
         t1 = tlist(i); t2 = tlist(i+1); t3 = tlist(i+2);
         x21 = mol.atxyz(1:3,t2) - mol.atxyz(1:3,t1);
         x23 = mol.atxyz(1:3,t2) - mol.atxyz(1:3,t3);
         pesc = dot(x21,x23);
         pvec = cross(x21,x23);
         pvec = sqrt(dot(pvec,pvec));
         angle = atan2(pvec,pesc)*180/pi;
         s1 = sprintf("%s(%d)", cell2mat(mol.atname(t1)), t1);
         s2 = sprintf("%s(%d)", cell2mat(mol.atname(t2)), t2);
         s3 = sprintf("%s(%d)", cell2mat(mol.atname(t3)), t3);
         printf("%12.6f %s-%s-%s\n", angle, s1, s2, s3);
      endfor
      printf("\n");
   endif
   if (tat > 3)
      printf("Dihedrals between quartets:\n");
      for i = 1:tat-3
         t1 = tlist(i); t2 = tlist(i+1); t3 = tlist(i+2); t4 = tlist(i+3);
         x12 = mol.atxyz(1:3,t2) - mol.atxyz(1:3,t1);
         x23 = mol.atxyz(1:3,t3) - mol.atxyz(1:3,t2);
         x34 = mol.atxyz(1:3,t4) - mol.atxyz(1:3,t3);
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
         s1 = sprintf("%s(%d)", cell2mat(mol.atname(t1)), t1);
         s2 = sprintf("%s(%d)", cell2mat(mol.atname(t2)), t2);
         s3 = sprintf("%s(%d)", cell2mat(mol.atname(t3)), t3);
         s4 = sprintf("%s(%d)", cell2mat(mol.atname(t4)), t4);
         printf("%12.6f %s-%s-%s-%s\n", dihedr, s1, s2, s3, s4);
      endfor
   endif

endfunction
