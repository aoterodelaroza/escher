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

function rep = rep_addcamera(repi,pos=-1,sky=[0 0 1],angle=45,persp=1)
% function rep = rep_addcamera(repi,pos,sky=[0 0 1],angle=45,persp=1)
%
% rep_addcamera - add a camera to a graphical representation
% with reasonable default values. pos is the camera location (by 
% default, oriented in the [0 0 1] direction and at a 
% distance to the object center) relative to the center of mass of the
% representation. angle controls the distance to the
% object (orthographic, persp=0) or the field of view (perspective,
% persp=1). sky is the sky vector (z by default).
%
% Input variables:
% repi: input representation.
% pos: position of the camera (3-element vector) relative to the
% representation center of mass.
% sky: direction of the sky vector.
% angle: distance to the object in terms of camera angle.
% persp: 1 for perspective, 0 for orthographic.
%
% Required output variables:
% rep: output representation.

  [xct xmin xmax xdel] = rep_getcm(repi);
  if (length(pos) != 3)
    pos = xct + [0 0 1] * xdel(3) * 2;
  else
    pos = xct + pos;
  endif

  ## set the camera 
  rep = repi;
  rep.cam = camera();
  rep.cam.location = pos; # position
  rep.cam.lookat = xct; # direction
  rep.cam.persp = persp; # perspective
  rep.cam.angle = angle; # field-of-view
  rep.cam.right = [1 0 0];
  rep.cam.up = [0 0 1];
  rep.cam.sky = sky;

endfunction
