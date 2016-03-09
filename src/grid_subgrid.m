% Copyright (c) 2015 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [newgrid] = grid_subgrid (grid, X0, X1, LOG=0)
% unction [newgrid] = grid_subgrid (grid, X0, X1, LOG=0)
%
% grid_subgrid - Returns a portion ofthe imput grid such that
%    all points in the new grid are contained in the [X0,X1]
%    parallelepiped.
%
% Required input variables:
% {grid}: Original grid structure.
% {X0,X1}: min and max values for the x, y, z, values of the pellelepiped
%   oriented by to the cartesian axes.
%
% Output:
% {newgrid}: new grid formed.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Cubefile structure according to gaussian:
% If the origin is (X0,Y0,Z0), and the increment is (X1,Y1,Z1), then point
% (I1,I2,I3) has the coordinates:
% 
%      X-coordinate: X0+(I1-1)*X1+(I2-1)*X2+(I3-1)*X3
%      Y-coordinate: Y0+(I1-1)*Y1+(I2-1)*Y2+(I3-1)*Y3
%      Z-coordinate: Z0+(I1-1)*Z1+(I2-1)*Z2+(I3-1)*Z3
%
% Althought this descriptions allows for a grid with non-orthogonal axis
% few or no implementation assumes this possibility. Therefore we assume the
% grid.dx matrix to be diagonal.
% Restriction: The new grid is a subset of the old one. Interpolation
% would be a more powerful method.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: Aug 2015
  #angtobohr = 1.88972613288564;
  #bohrtoans = 0.52917720859;

  #newgrid.f = grid.f( (grid.x0 + sum(grid.dx*diag([i,j,k])))>=X0 && (grid.x0 + sum(grid.dx*diag([i,j,k])))<=X1);

  xxg = zeros(size(grid.f));
  for i = 1:grid.n(1)
     for j = 1:grid.n(2)
        for k = 1:grid.n(3)
           xg = grid.x0 + sum(grid.dx*diag([i-1,j-1,k-1]));
           xxg(i,j,k) =  (xg'>=X0 && xg'<=X1);
        endfor
     endfor
  endfor

  newgrid.dx = grid.dx; # New and old grid share the interval sizes
  newgrid.f = grid.f(xxg==1);
  newgrid.n = size(newgrid.f);
  newgrid.omega = abs(det(newgrid.dx));
  #newgrid.a = diag(newgrid.n).*newgrid.dx;
  if (LOG>0)
     printf('grid_subgrid: Select a subset of a previous grid\n');
     printf('Old and new dimensions: %d %d %d and %d %d %d\n', grid.n, newgrid.n);
     printf('Old and new elem. num.: %dKb & %dKb\n', prod(size(grid.f))/1024, prod(size(newgrid.n))/1024);
  endif

endfunction
