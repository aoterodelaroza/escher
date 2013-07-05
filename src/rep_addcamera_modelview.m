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

function rep = rep_addcamera_modelview(repi,r,fov,LOG=0)
% rep = rep_addcamera_modelview(repi,r,fov,LOG=0)
%
% rep_addcamera_modelview - add a camera to a graphical representation, using
% the opengl modelview matrix.
% 
% Input variables:
% repi: input representation.
% r: opengl modelview matrix (4x4). Contains the model rotation (1:3,1:3)
%    and the translation (4,1:3). The (1:3,4) vector is the center of projection
%    (the point at which the camera points). If r is 4x3 instead of 4x4, then the 
%    cop is the barycenter of the balls in the representations.
% fov: opengl field of view (angle in pov).
% {LOG}: verbose level (0=silent,1=verbose).
%
% Required output variables:
% rep: output representation.

  ## xct = rep_getcm(repi);

  ## set the camera using euler angles
  rep = repi;
  rep.cam = struct();
  if (size(r,2) > 3)
    rep.cam.cop = r(1:3,4)'; ## use the cop in the modelview matrix
  else
    ## calculate the barycenter of the balls
    cm = [0 0 0];
    for i = 1:repi.nball
      cm += repi.ball{i}.x;
    endfor
    rep.cam.cop = cm / repi.nball;
  endif
  rep.cam.sky = [0 1 0]; ## default
  rep.cam.vuv = [0 1 0]; ## default
  rep.cam.rht = [1 0 0]; ## default divided by 4/3 (1:1 ratio)
  rep.cam.drt = [0 0 -1]; ## pov is right-handed
  rep.cam.angle = fov; ## field-of-view is angle
  rep.cam.matrix = r(1:4,1:3); ## modelview matrix
  rep.cam.invmatrix = 1; ## inverse of that

endfunction
