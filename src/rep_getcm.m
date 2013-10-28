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

function [xct xmin xmax xdel] = rep_getcm(rep)
% function rep_getcm(rep)
%
% rep_getcm - calculate the minimum, maximum, and center of mass of the
%             representation.
%
% Required input variables:
% rep: input representation.
%
% Output variables:
% xct: center of mass.
% xmin: minimum coordinates.
% xmax: maximum coordinates.
% xdel: xmax-xmin.

  ## calculate molecular limits
  xmax = zeros(1,3) - Inf;
  xmin = zeros(1,3) + Inf;
  for i = 1:rep.nball
    xmax = max(xmax,rep.ball{i}.x);
    xmin = min(xmin,rep.ball{i}.x);
  endfor
  for i = 1:rep.nstick
    xmax = max(xmax,rep.stick{i}.x0);
    xmin = min(xmin,rep.stick{i}.x0);
    xmax = max(xmax,rep.stick{i}.x1);
    xmin = min(xmin,rep.stick{i}.x1);
  endfor
  for i = 1:rep.nvertex
    xmax = max(xmax,rep.vertex{i}.x);
    xmin = min(xmin,rep.vertex{i}.x);
  endfor
  if (any(isinf(xmax)) || any(isinf(xmin)))
    xmin = xmax = xct = xdel = [0 0 0];
  else
    xct = 0.5 * (xmax+xmin);
    xdel = xmax - xmin + 1d-15;
  endif

endfunction
