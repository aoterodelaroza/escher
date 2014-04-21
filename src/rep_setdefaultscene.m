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

function rep = rep_setdefaultscene(repi,pos,sky=[0 0 1],lighting="direct",ifac=1,angle=45,persp=1);
% function rep = rep_setdefaultscene(repi,pos,sky=[0 0 1],lighting="direct",ifac=1,angle=45,persp=1);
%
% rep_setdefaultscene - given a representation, set up the camera
% using the camera position relative to the center of mass of the 
% representation. The sky vector can be set but is optional.
%
% Input variables:
% repi: input representation.
% pos: position of the camera relative to the representation center of mass.
% sky: direction of the sky vector.
% lighting: the lighting model.
% ifac: global light intensity factor
% angle: the distance from the object in terms of camera angle.
% persp: 1 for perspective, 0 for orthographic.
%
% Output variables:
% rep: output representation.

  rep = repi;

  ## add the camera, default placement
  rep = rep_addcamera(rep,pos,sky,angle,persp);

  ## add the lights, default placement
  rep = rep_lightmodel(rep,lighting,ifac);

  ## set the background color to white
  rep = rep_setbgcolor(rep,[255 255 255]);

endfunction
