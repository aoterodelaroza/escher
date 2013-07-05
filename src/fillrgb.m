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

function rgb = fillrgb(rgb0)
% function rgb = fillrgb(rgb0)
%
% fillrgb - fill a rgb triplet or quartet to a quintuplet using zeros.
%
% Required input variables:
% rgb0: input rgb triplet or quadruplet.
%
% Output:
% rgb: rgb quintuplet, same as input but filled with zeros.
%
  if (length(rgb0) == 3)
    rgb = [rgb0 0 0];
  elseif (length(rgb0) == 4)
    rgb = [rgb0 0];
  else
    rgb = rgb0;
  endif
endfunction
