% Copyright (c) 2012 Alberto Otero-de-la-Roza and Victor Lua~na
% Adapted from an AOR (2011) routine.
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

function c = cr_x2c(cr,x)
% function c = cr_x2c(cr,x)
%
% cr_x2c - crystallographic to cartesian coordinates. Use cr.r if available.
%          Otherwise, use the cholesky decomposition of the metric tensor.
%
% Required input variables:
% {cr}: crystal (structure)**2.
% {x}: crystallographic coordinates
%
% Output:
% {c}: cartesian coordinates.
%

  ## crystal to cartesian
  if (isfield(cr,"r"))
    r = cr.r;
  else
    if (isfield(cr,"g"))
      g = cr.g;
    else
      cc = cos(cr.b);
      g = cr.a' * cr.a;
      g(1,2) = g(2,1) = g(1,2) * cc(3);
      g(1,3) = g(3,1) = g(1,3) * cc(2);
      g(2,3) = g(3,2) = g(2,3) * cc(1);
    endif
    r = chol(g)';
  endif

  ## works for row and column
  if (size(x,2) > 1)
    c = x * r;
  else
    c = r' * x;
  endif

endfunction
