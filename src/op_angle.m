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

function a = op_angle(x0,x1,x2)
% function a = op_angle(x0,x1,x2)
%
% op_angle - returns the angle (degrees) between x0-x1-x2.

  n10 = x0-x1; n10 = n10/norm(n10);
  n12 = x2-x1; n12 = n12/norm(n12);
  a = acos(dot(n10,n12)) * 180/pi;

endfunction
