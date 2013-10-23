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

function rep = rep_setdefaultscene(repi,r="",angle=45,persp=1,LOG=0);
% function rep = rep_setdefaultscene(repi,r="",angle=45,persp=1,LOG=0);
%
% rep_setdefaultscene - given a representation, set up a camera, lights
% and background colors with reasonable default parameters.
%
% Input variables:
% repi: input representation.
% r: optional modelview matrix (4x3 or 4x4) or camangle vector (1x3).
% angle: the camera field of view.
% persp: 1 for perspective, 0 for orthographic.
% LOG: verbose level (0=silent,1=verbose).
%
% Output variables:
% rep: output representation.

  rep = repi;

  if (LOG > 0)
    printf("Adopting default parameters for representation...\n");
  endif

  ## add the camera, default placement
  if (isnumeric(r) && !isscalar(r))
    if (size(r,1) > 1)
      rep = rep_addcamera_modelview(rep,r,angle,persp,LOG);
    else
      rep = rep_addcamera(rep,r,angle,persp,LOG);
      r = "";
    endif
  else
    rep = rep_addcamera(rep,:,:,:,LOG);
    r = "";
  endif

  ## add the lights, default placement
  rep = rep_addlight(rep,rep.cam.cop,r,:,:,LOG);
  rep = rep_addlight(rep,rep.cam.sky*100,r,:,:,LOG);

  ## set the background color to white
  rep = rep_setbgcolor(rep,[255 255 255],LOG);

  if (LOG > 0)
    printf("\n");
  endif

endfunction
