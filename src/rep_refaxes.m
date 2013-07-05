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

function rep = rep_refaxes(addto="",scale=1,rad=0.01)
% function rep = rep_refaxes(addto="",scale=1,rad=0.01)
%
% rep_refaxes - add reference axes using 3 sticks (x=red, y=green, z=blue)
%               with length scale.
%
% Required input variables:
% addto: input representation.
% scale: length of the axes.
% rad: radius of the axes.
%
% Required output variables:
% rep: output representation.

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
  endif
  if (!isfield(rep,"nstick"))
    rep.nstick = 0;
    rep.stick = cell();
  endif

  ## center of mass
  [xct xmin xmax xdel] = rep_getcm(rep);

  ## add sticks
  n = rep.nstick;

  ## x axis
  n++;
  rep.stick{n}.name = "x-axis";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct + [1 0 0] * scale;
  rep.stick{n}.r = rad;
  rep.stick{n}.rgb = [255 0 0 0 0];
  rep.stick{n}.tex = "stick_default";
  n++;
  rep.stick{n}.name = "x-axis (negative)";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct - [1 0 0] * scale;
  rep.stick{n}.r = rad;
  rep.stick{n}.rgb = [28 0 0 0 0];
  rep.stick{n}.tex = "stick_default";
  ## y axis
  n++;
  rep.stick{n}.name = "y-axis";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct + [0 1 0] * scale;
  rep.stick{n}.r = rad;
  rep.stick{n}.rgb = [0 255 0 0 0];
  rep.stick{n}.tex = "stick_default";
  n++;
  rep.stick{n}.name = "y-axis (negative)";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct - [0 1 0] * scale;
  rep.stick{n}.r = rad;
  rep.stick{n}.rgb = [0 28 0 0 0];
  rep.stick{n}.tex = "stick_default";
  ## z axis
  n++;
  rep.stick{n}.name = "z-axis";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct + [0 0 1] * scale;
  rep.stick{n}.r = rad;
  rep.stick{n}.rgb = [0 0 255 0 0];
  rep.stick{n}.tex = "stick_default";
  n++;
  rep.stick{n}.name = "z-axis (negative)";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct - [0 0 1] * scale;
  rep.stick{n}.r = rad;
  rep.stick{n}.rgb = [0 0 28 0 0];
  rep.stick{n}.tex = "stick_default";
  rep.nstick = n;

endfunction
