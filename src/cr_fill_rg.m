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

function cr = cr_fill_rg(cr)
% function cr = cr_fill_rg(cr)
%
% cr_fill_rg - fill the r and g matrices using the cell parameters.
%
% Input:
% cr - the input crystal. This routine needs at least a and b. 
%
% Output:
% cr - 
%

  ## crystal to cartesian
  if (isfield(cr,"r"))
    if (!isfield(cr,"g"))
      cr.g = cr.r * cr.r';
    endif
  else
    if (!isfield(cr,"g"))
      cc = cos(cr.b);
      g = cr.a' * cr.a;
      g(1,2) = g(2,1) = g(1,2) * cc(3);
      g(1,3) = g(3,1) = g(1,3) * cc(2);
      g(2,3) = g(3,2) = g(2,3) * cc(1);
      cr.g = g;
    endif
    cr.r = chol(g)';
  endif

endfunction
