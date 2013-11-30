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

function x = cr_c2x(cr,c)
% function x = cr_c2x(cr,c)
%
% cr_c2x - cartesian to crystallographic coordinates. Use cr.r.
%
% Required input variables:
% {cr}: crystal (structure)**2.
% {c}: cartesian coordinates
%
% Output:
% {x}: crystallographic coordinates.
%

  ## crystal to cartesian
  rinv = inv(cr.r);

  ## works for row and column
  if (size(cz,2) > 1)
    x = c * rinv;
  else
    x = rinv' * c;
  endif

endfunction
