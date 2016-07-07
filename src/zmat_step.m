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

function xf = zmat_step(x0,x1,x2,d,ang,dieh)
% function xf = zmat_step(x0,x1,x2,d,ang,dieh)
%
% zmat_step - calculate the coordinates of atom number 4 given the coordinates
%             of atoms 1, 2 and 3, and the distance 4-1, angle 4-1-2 and dihedral
%             4-1-2-3
%
% Required input variables:
% x0..2: coordinates of the atoms.
% d: distance
% ang: angle (degrees)
% dieh: dihedral (degrees)
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: March 2012
  
  ## transform units
  ang *= pi / 180;
  dieh *= pi / 180;

  ## x0 is the origin
  x1 = x1 - x0;
  x2 = x2 - x0;

  ## x1-x0 is aligned to z
  crot(3,:) = x1 / norm(x1);
  crot(1,:) = cross(x1,[0 0 1]);
  if (abs(norm(crot(1,:))) < 1d-12)
    crot(1,:) = cross(x1,[0 1 0]);
  endif
  crot(1,:) = crot(1,:) / norm(crot(1,:));
  crot(2,:) = cross(crot(3,:),crot(1,:));
  
  ## transform x2
  x2p = crot * x2';
  [th2 ph2 r2] = cart2sph(x2p(1),x2p(2),x2p(3));
  rf = d;
  phf = pi/2 - ang;
  thf = th2 + dieh;
  [xf(1) xf(2) xf(3)] = sph2cart(thf,phf,rf);
  xf = inv(crot) * xf' + x0';

endfunction
