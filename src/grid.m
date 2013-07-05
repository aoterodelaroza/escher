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

function g = grid()
% function g = grid()
%
% grid - create an empty grid structure.
%
% Output:
% {g}: the empty grid structure with all the fields defined.
%

  grid.name = "";
  g.x0 = [0 0 0];
  g.dx = zeros(3,3);
  g.a = zeros(3,3);
  g.n = [0 0 0];
  g.f = [];
  g.omega = 0;

endfunction
