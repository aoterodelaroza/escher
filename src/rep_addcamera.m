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

function rep = rep_addcamera(repi,camangle=[80 75 45],zoom=3,LOG=0)
% function rep = rep_addcamera(repi,camangle=[80 75 45],zoom=3,LOG=0)
%
% rep_addcamera - add a camera to a graphical representation
%
% Required input variables:
% repi: input representation.
%
% Optional input variables (all have default values):
% {LOG}: verbose level (0=silent,1=verbose).
%
% Required output variables:
% rep: output representation.

  [xct xmin xmax xdel] = rep_getcm(repi);
  if (LOG>0)
    printf("Setting the camera using camangle: %.1f %.1f %.1f\n",camangle)
    printf("Zoom : %.1f\n",zoom);
    printf("Figure center (ang): %.10f %.10f %.10f \n",xct);
    printf("Extension (ang): %.10f %.10f %.10f \n",xdel);
  endif

  ## set the camera using camangle
  rep = repi;
  rep.cam = struct();
  camangle = camangle * pi / 180;
  visdis = zoom * max(xdel) / tan(camangle(3)/2);
  ss = sin(camangle(1));
  rep.cam.cop = xct + visdis * [ss*cos(camangle(2)), ss*sin(camangle(2)), cos(camangle(1))];
  rep.cam.sky = [0 0 1];
  rep.cam.vuv = [0 0 1];
  rep.cam.rht = [1 0 0];
  rep.cam.drt = [0 1 0] * zoom;
  rep.cam.vrp = xct;

endfunction
