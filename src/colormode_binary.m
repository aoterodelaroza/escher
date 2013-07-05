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

function rgb = colormode_binary(u,f,cp=[255 0 0],cm=[0 255 0])
% function rgb = colormode_binary(u,f,cp,cm)
%
% colormode_binary - color function generator for rep_surface. Given an 
%   array of (u,v) pairs (u) and the surface scalar function (f), 
%   returns the corresponding array of colors (rgb). The cp color is
%   assigned if f(u) >= 0 and the cm color if f(u) < 0.
%
% Required input variables:
% u: array (n,2) of points on the surface.
% f: z = f(u) scalar function.
% cp: positive color.
% cm: negative color.
%
% Required output variables:
% rgb: array (n,5) of colors at points u.
%

  rgb = zeros(size(u,1),5);
  cp = fillrgb(cp); cm = fillrgb(cm);
  for i = 1:size(u,1)
    if (f(u(i,:)) >= 0)
      rgb(i,:) = cp;
    else
      rgb(i,:) = cm;
    endif
  endfor

endfunction
