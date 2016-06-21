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

function d = op_dihedral(x0,x1,x2,x3)
% function d = op_dihedral(x0,x1,x2,x3)
%
% op_angle - returns the dihedral angle (degrees) between x0-x1-x2-x3.

  b1 = x1-x0; b2 = x2-x1; b3 = x3-x2;
  b12 = cross(b1, b2); b23 = cross(b2, b3);
  d = -atan2(norm(b2) * dot(b1, b23), dot(b12, b23)) * (180/pi);

endfunction
