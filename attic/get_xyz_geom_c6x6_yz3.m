#! /usr/bin/octave -q
% Copyright (C) 2010 Victor Lua~na and Alberto Otero-de-la-Roza
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

function geom = get_xyz_geom_c6x6_yz3 (filexyz,LOG=1)
% function geom = get_xyz_geom_c6x6_yz3 (filexyz,LOG=1)
%
% get_xyz_geom_c6x6_yz3 - Produces the significative geometry for the
% C6X6--YZ3 adduct.
%
% Required input variables:
% filexyz: name of the input data file in xyz format.
%
% % Optional input variables (all have default values):
% {LOG}: current level of printing. 0: no printing.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: November 2011

err = mol_dbstart(0);

mol = mol_readxyz(filexyz,:,0);

# Check the internal geometry?
mol = mol_internalgeometry(mol);

if (mol.atnumber(5) == 6)

   # Atoms: 1: Y; 2-4: Z_3; 5-10: C_6; 11-16: X_6
   # Center of mass of C_6 (column vector):
   cm = mean(mol.atxyz(1:3,5:10),2);

   # Normal vector for C_6 best fit plane:
   ###A = [mol.atxyz(1,5:10)'-cm(1), mol.atxyz(2,5:10)'-cm(2), mol.atxyz(3,5:10)'-cm(3)];
   A = (mol.atxyz(1:3,5:10)-cm(1:3,1)*[1,1,1,1,1,1])';
   [u,s,v] = svd(A);
   nc6 = v(:,end);
   # normalize?
   nc6 = nc6 / sqrt(dot(nc6,nc6));

   # Y-cm and Y-C_6(plane) distances:
   xx = mol.atxyz(1:3,1) - cm;
   dycm = sqrt(dot(xx,xx));
   %%%%dyc6 = dot(nc6, mol.atxyz(1:3,1));
   dyc6 = dot(nc6, xx);

   # Lateral displacement of the projection of Y with respect to the cm:
   l = sqrt(dycm*dycm - dyc6*dyc6);

   # Normal vector for the Z_3 group:
   cmZ3 = mean(mol.atxyz(1:3,2:4),2);
   A = (mol.atxyz(1:3,2:4)-cmZ3(1:3,1)*[1,1,1])';
   [u,s,v] = svd(A);
   nz3 = v(:,end);
   # normalize?
   nz3 = nz3 / sqrt(dot(nz3,nz3));

   # Angle between the nc6 and nz3 normals:
   angle = acosd(dot(nc6,nz3));

   # avg(YZ) dist and avg(ZYZ) angle
   ii = 0
   for i = 2:4
      yz(++ii) = mol_dist(mol, 1, i);
   endfor
   a_yz = mean(yz); s_yz = std(yz);
   ii = 0;
   zyz(++ii) = mol_angle(mol, 2, 1, 3);
   zyz(++ii) = mol_angle(mol, 3, 1, 4);
   zyz(++ii) = mol_angle(mol, 4, 1, 2);
   a_zyz = mean(zyz); s_zyz = std(zyz);

   # avg(CC), avg(CY), avg(CCC)
   ii = 0;
   t = 5 : 10;
   for i1 = 1 : 6
      i2 = mod(i1,6) + 1;
      i3 = mod(i2,6) + 1;
      cc(++ii) = mol_dist(mol, t(i1), t(i2));
      ccc(ii) = mol_angle(mol, t(i1), t(i2), t(i3));
      cx(ii) = mol_dist(mol, t(i1), t(i1)+6);
   endfor
   a_cc  = mean(cc);  s_cc  = std(cc);
   a_ccc = mean(ccc); s_ccc = std(ccc);
   a_cx  = mean(cx);  s_cx  = std(cx);

elseif (mol.atnumber(1) == 6)

   # Atoms: 1-6: C_6; 7-12; X_6; 13: Y; 14-16: Z_3.
   # Center of mass of C_6 (column vector):
   cm = mean(mol.atxyz(1:3,1:6),2);

   # Normal vector for C_6 best fit plane:
   A = (mol.atxyz(1:3,1:6)-cm(1:3,1)*[1,1,1,1,1,1])';
   [u,s,v] = svd(A);
   nc6 = v(:,end);
   # normalize?
   nc6 = nc6 / sqrt(dot(nc6,nc6));

   # Y-cm and Y-C_6(plane) distances:
   xx = mol.atxyz(1:3,13) - cm;
   dycm = sqrt(dot(xx,xx));
   dyc6 = dot(nc6, xx);

   # Lateral displacement of the projection of Y with respect to the cm:
   l = sqrt(dycm*dycm - dyc6*dyc6);

   # Normal vector for the Z_3 group:
   cmZ3 = mean(mol.atxyz(1:3,14:16),2);
   A = (mol.atxyz(1:3,14:16)-cmZ3(1:3,1)*[1,1,1])';
   [u,s,v] = svd(A);
   nz3 = v(:,end);
   # normalize?
   nz3 = nz3 / sqrt(dot(nz3,nz3));

   # Angle between the nc6 and nz3 normals:
   angle = acosd(dot(nc6,nz3));

   # avg(YZ) dist and avg(ZYZ) angle
   ii = 0
   for i = 14 : 16
      yz(++ii) = mol_dist(mol, 13, i);
   endfor
   a_yz = mean(yz); s_yz = std(yz);
   ii = 0;
   zyz(++ii) = mol_angle(mol, 14, 13, 15);
   zyz(++ii) = mol_angle(mol, 15, 13, 16);
   zyz(++ii) = mol_angle(mol, 16, 13, 14);
a_zyz = mean(zyz); s_zyz = std(zyz);

   # avg(CC), avg(CY), avg(CCC)
   ii = 0;
   for i1 = 1 : 6
      i2 = mod(i1,6) + 1;
      i3 = mod(i2,6) + 1;
      cc(++ii) = mol_dist(mol, i1, i2);
      ccc(ii) = mol_angle(mol, i1, i2, i3);
      cx(ii) = mol_dist(mol, i1, i1+6);
   endfor
   a_cc  = mean(cc);  s_cc  = std(cc);
   a_ccc = mean(ccc); s_ccc = std(ccc);
   a_cx  = mean(cx);  s_cx  = std(cx);

endif

geom.Ydist = dycm;
geom.Yproj = dyc6;
geom.d_lat = l;
geom.C6_Z3_angle = angle;
geom.YZ = a_yz;
geom.ZYZ = a_zyz;
geom.CC = a_cc;
geom.CX = a_cx;
geom.CCC = a_ccc;

if (LOG > 0)
   printf("Distance Y-c.mass(C_6) ........ %.5f\n", dycm);
   printf("Projection Y-on-C_6 ........... %.5f\n", dyc6);
   printf("Lateral displacement proj-cm .. %.5f\n", l);
   printf("C6 and Z3 planes angle (degs.). %.5f  %.5f\n", angle, 180-angle);
   printf("\n--Y_cm------l------YZ--------ZYZ--------CC------CX------CCC--\n");
   printf("%7.5f  %6.4f   %6.4f   %7.3f ", dycm, l, a_yz, a_zyz);
   printf(" %7.4f  %7.4f  %7.3f\n", a_cc, a_cx, a_ccc);
   printf("---------------- %7.4f  %8.4f ", s_yz, s_zyz);
   printf(" %7.4f  %7.4f %8.3f\n", s_cc, s_cx, s_ccc);
endif

endfunction

## Use it as a script
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      for i = 1 : nargin
         get_xyz_geom_c6x6_yz3(args{i});
      endfor
   else
      printf('Use as: get_xyz_geom_c6x6_yz3 file(s)\n');
   endif
endif
