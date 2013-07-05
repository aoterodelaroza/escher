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

function m = op_rot3D (a1,a2,a3, t = [0,0,0]', mode="euler_yxz", LOG=0)
% function m = op_rot3D (a1,a2,a3, t = [0,0,0]', mode="euler_yxz", LOG=0)
%
% op_rot3D - returns the matrix corresponding to a rotation in 3D defined
% in terms of three angles plus the translation contained in t.
%
% Required input variables:
% a1,a2,a3: the three rotation angles in degrees.
%
% Optional input variables (all have default values):
% {t = [0,0,0]'}: optional traslation to be done after the rotation.
% {mode = "euler_yxz"}: convention used to define the three rotation axes.
%       It can be any of:
%       * "euler_ABC" where ABC is a permutation of "xyz". It corresponds to
%         a sucession of counter-clockwise rotations: first rotA(a1), then
%         rotB(a2), and finally rotC(a3).
%       * "airplane" or "ypw": a1 is the yaw angle (gui~nada, i.e. rotation
%         around the vertical axis); a2 is the pitch angle (cabeceo, i.e.
%         rotation around the transverse axis, parallel to the wings of the
%         airplane); and a3 is the roll angle (alabeo, i.e. rotation around
%         the longitudinal or tail-to-nose axis of the airplane).
%       * "cwlh": clockwise/left-handed (phi,theta,psi) rotations.
%       * "axis": perform a rotation of a2 degrees around the unit axis
%         defined in the vector a1.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

if (strcmpi(mode,"euler_xyz"))
   m = op_rotz(a3) * op_roty(a2) * op_rotx(a1);
elseif (strcmpi(mode,"euler_xzy"))
   m = op_roty(a3) * op_rotz(a2) * op_rotx(a1);
elseif (strcmpi(mode,"euler_yxz"))
   m = op_rotz(a3) * op_rotx(a2) * op_roty(a1);
elseif (strcmpi(mode,"euler_yzx"))
   m = op_rotx(a3) * op_rotz(a2) * op_roty(a1);
elseif (strcmpi(mode,"euler_zxy"))
   m = op_roty(a3) * op_rotx(a2) * op_rotz(a1);
elseif (strcmpi(mode,"euler_zyx"))
   m = op_rotx(a3) * op_roty(a2) * op_rotz(a1);
elseif (strcmpi(mode,"airplane") | strcmpi(mode,"ypw"))
   m = op_rotx(a3) * op_roty(a2) * op_rotz(a1);
elseif (strcmpi(mode,"cwlh"))
   c1 = cos(a1*pi/180); s1 = sin(a1*pi/180);
   c2 = cos(a2*pi/180); s1 = sin(a2*pi/180);
   c3 = cos(a3*pi/180); s1 = sin(a3*pi/180);
   m = [ c2*c3, -c1*s3+s1*s2*c3,  s1*s3+c1*c2*c3; \
         c2*s3,  c1*c3+s1*s2*s3, -s1*c3+c1*s2*s3; \
        -s2,     s1*c2,           c1*c2           ];
elseif (strcmpi(mode,"axis"))
   # We check the normalization of the a2 vector and its row/column format:
   # p is the normalized column vector defining the rotation axis.
   st = sin(a2*pi/180); ct = cos(a2*pi/180); ct1 = 1 - ct;
   if (prod(size(a1)==[3,1]) == 1)
      p = a1 / (a1'*a1);
   elseif (prod(size(a1)==[1,3]) == 1)
      p = a1' / (a1*a1');
   else
      error('rot_3D/axis: a1 should be a 3x1 unit vector!',
   endif
   m = (p * p') * ct1 + \
     [ +ct,      +p(3)*st, -p(2)*st; \
       -p(3)*st, +ct,      +p(1)*st; \
       +p(2)*st, -p(1)*st, +ct       ];
else
   error('op_rot3D: Unknown 3 angles rotation mode!');
endif

[r,c] = size(t);
if (r==3 & c==1)
   m(:,4) = t;
elseif (r==1 & c==3)
   m(:,4) = t';
else
   error('op_rot3D: wrong translation component!');
endif

if (LOG > 0)
   printf("op_rot3D: %s rotation\n", mode);
   if (strcmpi(mode,"axis"))
      printf("Rotation axis vector:", %.3f %.3f %.3f\n", a1);
      printf("Rotation angle (deg):", %.3f\n", a2);
   else
      printf("Angles (deg) :", %.3f %.3f %.3f\n", a1, a2, a3);
   endif
   printf("Traslation vector   :", %.3f %.3f %.3f\n", t);
   printf("Rotation matrix (Op * column_vector):\n");
   printf("%12.5e %12.5e %12.5e %12.5e\n", m);
endif

endfunction
