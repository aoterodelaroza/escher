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

function rep = rep_refaxes(addto="",refrep="",scale=1,srad=1)
% function rep = rep_refaxes(addto="",refrep="",scale=1,srad=1)
%
% rep_refaxes - add reference axes using 3 sticks (x=red, y=green, z=blue)
%               with length scale.
%
% Required input variables:
% addto: input representation.
% refrep: representation that gives the center of mass for the axes.
% scale: length of the axes.
% srad: radius of the axes is srad times 0.05.
%
% Required output variables:
% rep: output representation.

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation_();
  endif

  ## center of mass
  if (!isempty(refrep) && isstruct(refrep)) 
      [xct xmin xmax xdel] = rep_getcm(refrep);
  else
      xct = [0 0 0];
  endif

  ## add sticks
  n = rep.nstick;

  ## register textures in the representations
  [rep itex] = rep_registertexture(rep,"stick_default");

  ## x axis
  n++;
  rep.stick{n} = stick();
  rep.stick{n}.name = "x-axis";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct + [1 0 0] * scale;
  rep.stick{n}.r = srad * 0.05;
  rep.stick{n}.rgb = [255 0 0 0 0];
  rep.stick{n}.tex = itex;

  ## y axis
  n++;
  rep.stick{n} = stick();
  rep.stick{n}.name = "y-axis";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct + [0 1 0] * scale;
  rep.stick{n}.r = srad * 0.05;
  rep.stick{n}.rgb = [0 255 0 0 0];
  rep.stick{n}.tex = itex;

  ## z axis
  n++;
  rep.stick{n} = stick();
  rep.stick{n}.name = "z-axis";
  rep.stick{n}.x0 = xct + [0 0 0];
  rep.stick{n}.x1 = xct + [0 0 1] * scale;
  rep.stick{n}.r = srad * 0.05;
  rep.stick{n}.rgb = [0 0 255 0 0];
  rep.stick{n}.tex = itex;

  rep.nstick = n;

endfunction
