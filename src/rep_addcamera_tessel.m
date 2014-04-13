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

function rep = rep_addcamera_tessel(repi,camangle=[80 75 45],angle=45,persp=1,lookat=-1)
% function rep = rep_addcamera_tessel(repi,camangle=[80 75 45],angle=45,persp=1,lookat=-1)
%
% rep_addcamera_tessel - add a camera to a graphical representation
%  using tessel's positioning method.
%
% Input variables:
% repi: input representation.
% camangle: the three camera angles in spherical coordinate.s The first is
%           phi (angular separation from the +z axis), the second is theta
%           (rotation about the z) and the third controls the distance
%           from the center of mass. 
% angle: distance to the object in terms of camera angle.
% persp: 1 for perspective, 0 for orthographic.
% lookat: 3d vector specifying the point at which the camera is pointed.
%         By default (-1), it is the barycenter of the representation.
%
% Output variables:
% rep: output representation.

  [xct xmin xmax xdel] = rep_getcm(repi);
  if (length(lookat) != 3)
    lookat = xct;
  endif

   ## set the camera using camangle
   rep = repi;
   rep.cam = camera();
 
   camangle = camangle * pi / 180;
   visdis = max(xdel) / tan(camangle(3)/2);
   ss = sin(camangle(1));
   rep.cam.location = xct + visdis * [ss*cos(camangle(2)), ss*sin(camangle(2)), cos(camangle(1))];
   rep.cam.angle = angle;
   rep.cam.persp = persp;
   rep.cam.lookat = xct;
   rep.cam.right = [1 0 0];
   rep.cam.up = [0 0 1];
   rep.cam.sky = [0 0 1];

endfunction
