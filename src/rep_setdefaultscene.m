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

function rep = rep_setdefaultscene(repi,r="",lighting="3point",ifac=1,angle=45,persp=1);
% function rep = rep_setdefaultscene(repi,r="",lighting="3point",ifac=1,angle=45,persp=1);
%
% rep_setdefaultscene - given a representation, set up a camera, lights
% and background colors with reasonable default parameters.
%
% Input variables:
% repi: input representation.
% r: optional camangle vector (1x3) for use with rep_addcamera_tessel.
% lighting: the lighting model.
% ifac: global light intensity factor
% angle: the distance from the object in terms of camera angle.
% persp: 1 for perspective, 0 for orthographic.
%
% Output variables:
% rep: output representation.

  rep = repi;

  ## add the camera, default placement
  if (isnumeric(r) && !isscalar(r))
    rep = rep_addcamera_tessel(rep,r,angle,persp);
  else
    rep = rep_addcamera(rep);
  endif

  ## add the lights, default placement
  rep = rep_lightmodel(rep,lighting,ifac);

  ## set the background color to white
  rep = rep_setbgcolor(rep,[255 255 255]);

endfunction
