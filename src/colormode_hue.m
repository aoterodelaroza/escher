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

function rgb = colormode_hue(u,f)
% function rgb = colormode_hue(u,f)
%
% colormode_hue - color function generator for rep_surface. Given an 
%   array of (u,v) pairs (u) and the surface scalar function (f), 
%   returns the corresponding array of colors (rgb). The color is calculated
%   using the value of f as hue in the hsv scheme. f should give values in the 
%   range 0 to 360, with {0,60,120,180,240,300} being interpreted as
%   {red, yellow, green, cyan, blue, magenta}.
%
% Required input variables:
% u: array (n,2) of points on the surface.
% f: z = f(u) scalar function.
%
% Required output variables:
% rgb: array (n,5) of colors at points u.
%

  fval = f(u);
  rgb = zeros(size(u,1),5);
  for i = 1:size(u,1)
    h1 = mod(floor(fval(i)),360) / 60;
    ii = floor(h1);
    ff = h1 - ii;
    s = 1; v = 1; p = 0; q = 1 - ff; t = ff;
    if (ii == 0)
      rgb(i,1:3) = [v t p];
    elseif (ii == 1)
      rgb(i,1:3) = [q v p];
    elseif (ii == 2)
      rgb(i,1:3) = [p v t];
    elseif (ii == 3)
      rgb(i,1:3) = [p q v];
    elseif (ii == 4)
      rgb(i,1:3) = [t p v];
    elseif (ii == 5)
      rgb(i,1:3) = [v p q];
    endif
  endfor

endfunction
