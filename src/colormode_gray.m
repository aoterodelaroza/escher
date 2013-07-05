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

function rgb = colormode_gray(u,f,scale=[])
% function rgb = colormode_gray(u,f)
%
% colormode_gray - color function generator for rep_surface. Given an 
%   array of (u,v) pairs (u) and the surface scalar function (f), 
%   returns the corresponding array of colors (rgb). The scale is calculated
%   by default using the minimum and maximum of f on the u list, unless
%   scale = [fmin fmax] is given. The color scheme is grayscale.
%
% Required input variables:
% u: array (n,2) of points on the surface.
% f: z = f(u) scalar function.
% scale: [fmin fmax] values for the grayscale.
%
% Required output variables:
% rgb: array (n,5) of colors at points u.
%

  fval = f(u);
  if (isempty(scale))
    scale = [min(fval), max(fval)];
  endif

  fdel = 2 / (scale(2)-scale(1));
  rgb = zeros(size(u,1),5);
  for i = 1:size(u,1)
    ff = (fval(i)-scale(1)) * fdel - 1;
    rgb(i,1:3) = 0.5 * (1 - ff);
  endfor

endfunction
