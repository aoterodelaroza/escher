#! /etc/alternatives/octave -q
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

function [cell,G,V] = cr_pw2cell(ibrav, celldm, cell_param, LOG=1)
%function cell = cr_pw2cell(ibrav, celldm, cell_param, LOG=1)
%
% cr_pw2cell - convert the crystal output description produced by pwscf
% into standard crystallographic cell parameters.
%
% Required input variables:
% {ibrav}: the Bravais cell defined in pwscf.
% {celldm()}: the 6 components vector containin cell data.
% {cell_param(,)}: the three v1..v3 vectors row organized into a 6x6 matrix.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {cell}: crystallographic (a,b,c,alpha,beta,gamma) description of the cell.
%
% Optional output variables:
% {G}: metrical matrix.
% {V}: cell volume.
%
% Referencia: pwscf, routine flib/latgen.f90 and documentation of the code.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: April 2012

switch (ibrav)
  case {0,13,14,-5,-12}
    R = celldm(1) * cell_param;
    G = R*R';
    V = sqrt(abs(det(G)));
    a = sqrt(G(1,1));
    b = sqrt(G(2,2));
    c = sqrt(G(3,3));
    alpha = acos(R(2,3)/(b*c)) * 180/pi;
    beta  = acos(R(3,3)/(a*c)) * 180/pi;
    gamma = acos(R(1,2)/(a*b)) * 180/pi;
    cell = [a,b,c,alpha,beta,gamma];
    return;
  case 1
    # cubic P (cP, sc)
    R = celldm(1) * cell_param;
    G = R*R';
    V = sqrt(abs(det(G)));
    a = V**(1/3);
    cell = [a,a,a,90,90,90];
    ok = compare(R, a*eye(3));
    return;
  case 2
    # cubic F (cF, fcc)
    R = celldm(1) * cell_param;
    G = R*R';
    V = sqrt(abs(det(G)));
    a = (4*V)**(1/3);
    cell = [a,a,a,90,90,90];
    ok = compare(R, (a/2)*[-1,0,1; 0,1,1; -1,1,0]);
    return;
  case 3
    # cubic I (cI, bcc)
    R = celldm(1) * cell_param;
    G = R*R';
    V = sqrt(abs(det(G)));
    a = (4*V)**(1/3);
    cell = [a,a,a,90,90,90];
    ok = compare(R, (a/2)*[1,1,1; -1,1,1; -1,-1,1]);
    return;
  case 4
    # hexagonal and trigonal P (hP)
    R = celldm(1) * cell_param;
    G = R*R';
    V = abs(det(R));
    a = R(1,1);
    c = R(3,3);
    ok = compare(R, a*[1,0,0; -1/2,sqrt(3)/2,0; 0,0,c/a]);
    cell = [a,a,c,90,90,120];
    return;
  case 5
    # trigonal R (hR) c is the threefold axis
    R = celldm(1) * cell_param;
    G = R*R';
    rr = ( R(1,1) / R(1,3) )**2;
    cc = (3-2*rr) / ( 4*rr + 3);
    alpha = acos(cc) * 180/pi;
    a = R(1,1) / sqrt((1-cc)/3);
    tx = sqrt((1-cc)/2);
    ty = sqrt((1-cc)/6);
    tz = sqrt((1+2*cc)/3);
    ok = compare(R, a*[tx,-ty,tz; 0,2*ty,tz; -tx,-ty,tz]);
    cell = [a,a,a,alpha,alpha,alpha];
  case 6
    # tetragonal P (tP, st)
    R = celldm(1) * cell_param;
    G = R*R';
    a = R(1,1);
    c = R(3,3);
    ok = compare(R, [a,0,0; 0,a,0; 0,0,c]);
    cell = [a,a,c,90,90,90];
  case 7
    R = celldm(1) * cell_param;
    G = R*R';
    V = abs(det(R));
    a = 2 * R(1,1);
    c = 2 * R(1,3);
    ok = compare(R, (a/2) * [1,-1,c/a; 1,1,c/a; -1,-1,c/a]);
    cell = [a,a,c,90,90,90];
    return;
  case 8
    # orthorhombic P (oP, so)
    R = celldm(1) * cell_param;
    G = R*R';
    a = R(1,1);
    b = R(2,2);
    c = R(3,3);
    ok = compare(R, [a,0,0; 0,b,0; 0,0,c]);
    cell = [a,b,c,90,90,90];
  case 9
    # orthorhombic c base centered (oC, cbco)
    R = celldm(1) * cell_param;
    G = R*R';
    a = 2 * R(1,1);
    b = 2 * R(1,2);
    c = R(3,3);
    ok = compare(R, [a/2,b/2,0; -a/2,b/2,0; 0,0,c]);
    cell = [a,b,c,90,90,90];
  case 10
    # face centered orthorhombic (oF, fco)
    R = celldm(1) * cell_param;
    G = R*R';
    a = 2 * R(1,1);
    b = 2 * R(2,2);
    c = 2 * R(3,3);
    ok = compare(R, [a/2,0,c/2; a/2,b/2,0; 0,b/2,c/2]);
    cell = [a,b,c,90,90,90];
  case 11
    # body centered orthorhombic (oI, bco)
    R = celldm(1) * cell_param;
    G = R*R';
    a = 2 * R(1,1);
    b = 2 * R(2,2);
    c = 2 * R(3,3);
    ok = compare(R, [a,b,c; -a,b,c; -a,-b,c]/2);
    cell = [a,b,c,90,90,90];
  case 12
    # monoclinic P (mP, sm)
    R = celldm(1) * cell_param;
    G = R*R';
    a = R(1,1);
    b = sqrt(R(2,1)**2 + R(2,2)**2);
    c = R(3,3);
    gamma = atan2(R(2,2)/b, R(2,1)/b) * 180/pi;
    ok = compare(R, [a,0,0; b*cosd(gamma),b*sind(gamma),0; 0,0,c]);
    cell = [a,b,c,90,90,gamma];
  otherwise
    cell = [1,1,1,90,90,90];
    error ('cr_pw2cell: Sorry, not yet implemented!');
endswitch
endfunction

function ok = compare (R, newR)
   if (sum(sum(abs(R-newR))) > 1e-6)
      ok = false;
      printf("warning: initial lattice is not reproduced!\n");
      printf("Inversion of G would produce:\n");
      G = R*R';
      a = sqrt(G(1,1));
      b = sqrt(G(2,2));
      c = sqrt(G(3,3));
      alpha = acos(R(2,3)/(b*c)) * 180/pi;
      beta  = acos(R(3,3)/(a*c)) * 180/pi;
      gamma = acos(R(1,2)/(a*b)) * 180/pi;
      printf(" %.6f  %.6f  %.6f", a, b, c);
      printf(" %.6f  %.6f  %.6f", alpha, beta, gamma);
      printf("\n");
   else
      ok = true;
   endif
endfunction

function [ibrav,celldm,axes] = cr_read(file)
[fid,msg] = fopen(file,"r");
if (fid < 0 || ferror(fid)) 
   disp(msg)
   error("cr_pw2cell(cr_read): Could not find -- %s",file);
endif

line = fgetl(fid);
[g{1:1},cout] = sscanf(line,"%s",'C');
ibrav = cellfun('str2num',g(1));

celldm = zeros(1,6);
line = fgetl(fid);
[g{1:1},cout] = sscanf(line,"%s",'C');
celldm(1) = cellfun('str2num',g(1));

axes = zeros(3,3);
for i = 1:3
   line = fgetl(fid);
   [g{1:3},cout] = sscanf(line,"%s %s %s",'C');
   axes(i,1:3) = cellfun('str2num',g(1:3));
endfor

endfunction

##function [g] = leenumeros(fid,n)
##endfunction

## Use it as a script
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      [ibrav,celldm,axes] = cr_read(args{1});
      cell = cr_pw2cell(ibrav,celldm,axes)
   else
      printf('To bo used from pw2cell.awk as: pw2cell.awk pwoutput\n');
   endif
endif
