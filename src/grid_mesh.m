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
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2012

  if (abs(g.a(1,2))+abs(g.a(1,3))+abs(g.a(2,3)+abs(g.a(2,1))+abs(g.a(3,1))+abs(g.a(3,2))) < 1e-14)
    x = g.x0(1) + g.dx(1,1) * (0:g.n(1)-1);
    y = g.x0(2) + g.dx(2,2) * (0:g.n(2)-1);
    z = g.x0(3) + g.dx(3,3) * (0:g.n(3)-1);
    [xx,yy,zz] = ndgrid(x,y,z);
  else
    xx = yy = zz = zeros(size(g.f));
    for i1 = 0:g.n(1)-1
      ii1 = i1+1;
      for i2 = 0:g.n(2)-1
        ii2 = i2+1;
        for i3 = 0:g.n(3)-1
          ii3 = i3+1;
          xx(ii1,ii2,ii3) = g.x0(1) + i1 * g.dx(1,1) + i2 * g.dx(2,1) + i3 * g.dx(3,1);
          yy(ii1,ii2,ii3) = g.x0(2) + i1 * g.dx(1,2) + i2 * g.dx(2,2) + i3 * g.dx(3,2);
          zz(ii1,ii2,ii3) = g.x0(3) + i1 * g.dx(1,3) + i2 * g.dx(2,3) + i3 * g.dx(3,3);
        endfor
      endfor
    endfor
  endif

endfunction
