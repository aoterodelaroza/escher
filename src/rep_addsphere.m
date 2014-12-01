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

function rep  = rep_addsphere(addto="", x0=[0 0 0], rad=1, wire=1,...
                              frgb=[0 0 128 115 115],ftex="opaque_triangle_default")
% function rep  = rep_addsphere(addto="", x0=[0 0 0], rad=1, wire=1,...
%                               frgb=[0 0 128 115 115],ftex="opaque_triangle_default")
%
% rep_addsphere - add a sphere to the representation.
%
% Input variables:
% addto: representation where the cube will be added. If none is given, a new rep is created.
% x0: center of the sphere
% rad: radius of the sphere
% wire: use a wireframe.
% frgb: rgbft for the face color (from 0 to 255). Make it empty ("" or []) to 
%       deactivate faces.
% ftex: face texture.
%
% Output variables:
% rep: the representation containing the cube.
%

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
    if (isfield(molc,"name") && !isempty(molc.name))
      rep.name = molc.name;
    endif
  endif

  ## register the textures, if present
  [rep iftex] = rep_registertexture(rep,ftex);

  ## build the cube using the convex hull
  n = rep.nball = rep.nball + 1;
  rep.ball{n}.x = x0;
  rep.ball{n}.name = "user_ball";
  rep.ball{n}.r = rad;
  rep.ball{n}.rgb = fillrgb(frgb);
  rep.ball{n}.tex = iftex;
  rep.ball{n}.wire = wire;
  ## register shapes3.inc for loading if using rounded cylinders
  if (wire)
    rep.load.shapes3 = 1;
  endif

endfunction
