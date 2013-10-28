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

function rep = representation()
% function rep = representation()
%
% representation - create an empty rep structure.
%
% Output:
% {rep}: the empty representation.
%

  rep.name = "";

  rep.nball = 0;
  rep.ball = cell();

  rep.nstick = 0;
  rep.stick = cell();

  rep.ntriangle = 0;
  rep.nvertex = 0;
  rep.triangle = cell();
  rep.vertex = cell();

  rep.cam = camera();
  rep.nlight = 0;
  rep.light = cell();
  rep.bgcolor = zeros(1,3);

endfunction
