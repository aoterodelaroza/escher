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

function rep = rep_addcamera_view3dscene(repi,pos,dir,up,angle=45,persp=1)
% function rep = rep_addcamera_view3dscene(repi,pos,dir,up,angle=45,persp=1)
%
% rep_addcamera_view3dscene - add a camera to a graphical representation
% using the location, lookat, and sky vectors. These can be obtained
% from the most recent version of view3dscene. They appear written on
% the visualization window, under the names "pos", "dir", and "up",
% respectively. The obtained povray plot is in the same orientation as
% in view3dscene. angle controls the distance to the object
% (orthographic, persp=0) or the field of view (perspective, persp=1).
%
% Input variables:
% repi: input representation.
% pos: position of the camera (3-element vector).
% dir: direction of the camera (3-element vector).
% up:  up vector of the camera (3-element vector).
% angle: distance to the object in terms of camera angle.
% persp: 1 for perspective, 0 for orthographic.
%
% Required output variables:
% rep: output representation.

  ## set the camera 
  rep = repi;
  rep.cam = camera();
  rep.cam.location = pos; # position
  rep.cam.lookat = pos + 10 * dir; # direction
  rep.cam.persp = persp; # perspective
  rep.cam.angle = angle; # field-of-view
  rep.cam.right = [1 0 0];
  rep.cam.up = [0 0 1];
  rep.cam.sky = up;

endfunction
