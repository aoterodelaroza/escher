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

function rep  = rep_addcube(addto="", x0=[0 0 0], x1=[1 1 1], frgb=[0 0 128 115 115], ergb=[0 0 128],...
                            ftex="opaque_triangle_default", etex="stick_default",erad=0.025)
% function [rep molc1 molv1]  = rep_addcube(addto="", x0=[0 0 0], x1=[1 1 1], frgb=[0 0 128 115 115], ergb=[0 0 128],...
%                                           ftex="opaque_triangle_default", etex="stick_default",erad=0.025)
%
% rep_addcube - add a cube to the representation.
%
% Input variables:
% addto: representation where the cube will be added. If none is given, a new rep is created.
% x0, x1: limits of the parallelepiped
% frgb: rgbft for the face color (from 0 to 255). Make it empty ("" or []) to 
%       deactivate faces.
% ergb: rgbft for the edge color (from 0 to 255). Make it empty ("" or []) to 
%       deactivate edges. 
% ftex: face texture.
% etex: edge texture.
% erad: edge radius.
%
% Output variables:
% rep: the representation containing the cube.
%

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation_();
    if (isfield(molc,"name") && !isempty(molc.name))
      rep.name = molc.name;
    endif
  endif

  ## register the textures, if present
  if (!isempty(frgb))
    [rep iftex] = rep_registertexture(rep,ftex);
  endif
  if (!isempty(ergb))
    [rep ietex] = rep_registertexture(rep,etex);
  endif

  ## build the cube using the convex hull
  v = [x0(1) x0(2) x0(3)
       x1(1) x0(2) x0(3)
       x0(1) x1(2) x0(3)
       x0(1) x0(2) x1(3)
       x1(1) x1(2) x0(3)
       x1(1) x0(2) x1(3)
       x0(1) x1(2) x1(3)
       x1(1) x1(2) x1(3)
       ];
  h = convhulln(v);
  if (!isempty(frgb))
    nv0 = rep.nvertex;
    for i = 1:8
      rep.nvertex += 1;
      rep.vertex{rep.nvertex} = vertex();
      rep.vertex{rep.nvertex}.x = v(i,:);
      rep.vertex{rep.nvertex}.rgb = fillrgb(frgb);
    endfor
    for i = 1:size(h,1)
      rep.ntriangle += 1;
      rep.triangle{rep.ntriangle} = triangle();
      rep.triangle{rep.ntriangle}.idx = nv0 + h(i,:);
      rep.triangle{rep.ntriangle}.rgb = fillrgb(frgb);
      rep.triangle{rep.ntriangle}.tex = iftex;
    endfor
  endif
  if (!isempty(ergb))
    icon = zeros(8);
    kk = [1 2; 1 3; 2 3];
    for i = 1:size(h,1)
      for j = 1:3
        rep.nstick = rep.nstick + 1;
        rep.stick{rep.nstick} = stick();
        rep.stick{rep.nstick}.name = "";
        rep.stick{rep.nstick}.x0 = v(h(i,kk(j,1)),:);
        rep.stick{rep.nstick}.x1 = v(h(i,kk(j,2)),:);
        rep.stick{rep.nstick}.r = erad;
        rep.stick{rep.nstick}.rgb = ergb;
        rep.stick{rep.nstick}.tex = ietex;
      endfor
    endfor
  endif

endfunction
