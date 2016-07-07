% Copyright (c) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [xx, yy, zz] = grid_mesh(g)
% function [xx, yy, zz] = grid_mesh(g)
%
% grid_mesh - given a grid of function values, calculate xx, yy and zz
% grids of the same shape containing the x, y and z coordinates. Equivalent
% to ndgrid (an xy-transpose of meshgrid).
%
% Required input variables:
% g: the grid.
%
% Output:
% xx: a three-dimensional array the same shape as g.f containing the
% x coordinates.
% yy, zz: same as above.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: Jan 2012

  if (abs(g.a(1,2))+abs(g.a(1,3))+abs(g.a(2,3)+abs(g.a(2,1))+abs(g.a(3,1))+abs(g.a(3,2))) < 1e-14)
    x = g.x0(1) + g.dx(1,1) * (0:g.n(1)-1);
    y = g.x0(2) + g.dx(2,2) * (0:g.n(2)-1);
    z = g.x0(3) + g.dx(3,3) * (0:g.n(3)-1);
    [xx,yy,zz] = ndgrid(x,y,z);
  else
    [x,y,z] = ndgrid(1:g.n(1),1:g.n(2),1:g.n(3));
    x -= 1; y -= 1; z -= 1;
    xx = g.x0(1) + x * g.dx(1,1) + y * g.dx(2,1) + z * g.dx(3,1);
    yy = g.x0(1) + x * g.dx(1,2) + y * g.dx(2,2) + z * g.dx(3,2);
    zz = g.x0(1) + x * g.dx(1,3) + y * g.dx(2,3) + z * g.dx(3,3);
  endif

endfunction
