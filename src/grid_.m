% Copyright (C) 2011 Alberto Otero-de-la-Roza and Victor Lua~na
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

function g = grid_()
% function g = grid_()
%
% grid - create an empty 3D grid structure. This is the grid constructor.
%        The structure corresponds to the cube format definition in gaussian
%        and it could correspond, in principle, to non-orthogonal grids.
%
% Output:
% {g}: the empty grid structure with all the fields defined.
%
% Description of the fields:
  g.name = "";
  g.x0 = [0 0 0];    # Origin coordinates for the grid points
  g.dx = zeros(3,3); # interval sizes
  g.a = zeros(3,3);  # range (diag(g.n) .* g.dx)
  g.n = [0 0 0];     # Number of points in each of the 3D directions
  g.omega = 0;       # Volume of a voxel element: abs(det(g.dx))
  g.f = [];          # Scalar field values at the grid points

endfunction
