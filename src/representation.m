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

function rep = representation_()
% function rep = representation_()
%
% representation_ - Vonstructor: create an empty rep structure.
%
% Output:
% {rep}: the empty representation.
%

  rep.name = "";

  rep.nball = 0;      # number of balls
  rep.ball = cell();  # Read the ball rep in mol_ball.m

  rep.nstick = 0;     # number of sticks
  rep.stick = cell(); # Raed the stick rep in mol_stick.m

  rep.ntriangle = 0;  # Triangles in polyhedra (See mol_polyhedron.m)
  rep.nvertex = 0;
  rep.triangle = cell();
  rep.vertex = cell();

  rep.cam = struct();        # Definition of the camera
  rep.nlight = 0;            # Number of lights in the scene
  rep.light = cell();        # Definition for each light
  rep.bgcolor = zeros(1,3);a # color in the backgound

endfunction
