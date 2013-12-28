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

function [rep surf] = grid_isosurface(g,addto="",iso,frgb=[51 160 255 0 0],ergb=[],...
                               ftex="opaque_triangle_default",etex="stick_default",erad=0.005)
% function [rep surf] = grid_isosurface(g,addto="",iso,frgb=[51 160 255 0 0],ergb=[],...
%                                ftex="opaque_triangle_default",etex="stick_default",erad=0.005)
%
% grid_isosurface - plot an isosurface into a representation.
%
% Input variables:
% g: grid description
% iso: isovalue.
%
% Optional input variables (all have default values):
% addto: input representation on which the new representation is added.
% frgb: the faces color. Make it empty to deactivate faces.
% ergb: the edges color. Make it empty to deactivate edges.
% ftex: face texture.
% etex: edge texture.
% erad: edge radius.
%
% Output:
% rep: the reprsentation containing the isosurface.
% surf: the surface object containing the isosurface.
  
  bohrtoans = 0.52917720859;

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
    if (isfield(g,"name") && !isempty(g.name))
      rep.name = g.name;
    endif
  endif

  ## build the isosurface
  surf = surface_();
  [xx,yy,zz] = grid_mesh(g);
  [surf.f,surf.v] = isosurface(xx,yy,zz,g.f,iso);
  surf.v *= bohrtoans;
  nor = cross((surf.v(surf.f(:,2),:) - surf.v(surf.f(:,1),:)),(surf.v(surf.f(:,3),:) - surf.v(surf.f(:,1),:)));
  surf.n = zeros(size(surf.v));
  for i = 1:rows(surf.f)
    surf.n(surf.f(i,1),:) += nor(i,:);
    surf.n(surf.f(i,2),:) += nor(i,:);
    surf.n(surf.f(i,3),:) += nor(i,:);
  endfor
  surf.n = surf.n ./ norm(surf.n,2,"rows");
  surf.n(find(!isfinite(surf.n))) = 1/sqrt(3);

  ## face texture
  if (!isempty(frgb))
    surf.frgb = fillrgb(frgb);
    [rep iftex] = rep_registertexture(rep,ftex);
    surf.ftex = iftex;
  endif

  ## add it to the representation
  rep.nsurf += 1;
  rep.surf{rep.nsurf} = surf;

#  if (!isempty(frgb))
#    frgb = fillrgb(frgb);
#
#    ## pre-allocate vertices
#    nv = rep.nvertex;
#    rep.nvertex += size(v,1);
#    empty = struct("x",[0 0 0],"rgb",frgb);
#    rep.vertex(nv+1:rep.nvertex) = {empty};
#    ## fill
#    for i = nv+1:rep.nvertex
#      rep.vertex{i}.x = v(i-nv,:);
#    endfor
#
#    ## pre-allocate triangles
#    n = rep.ntriangle;
#    rep.ntriangle += size(f,1);
#    [rep iftex] = rep_registertexture(rep,ftex);
#    empty = struct("idx",[0 0 0],"rgb",[0 0 0 0 0],"tex",iftex);
#    rep.triangle(n+1:rep.ntriangle) = {empty};
#    ## fill triangles
#    for i = n+1:rep.ntriangle
#      rep.triangle{i}.idx = nv + f(i-n,:);
#      rep.triangle{i}.rgb = (rep.vertex{rep.triangle{i}.idx(1)}.rgb +...
#                             rep.vertex{rep.triangle{i}.idx(2)}.rgb +...
#                             rep.vertex{rep.triangle{i}.idx(3)}.rgb) /3;
#    endfor
#  endif

  ## edges
  if (!isempty(ergb))
    ## pre-determine the edge pairs 
    nt = size(surf.f,1);
    ipairs = zeros(3*nt,2);
    jaux = ((1:nt) - 1) * 3 + 1;
    for i = 1:nt
      aux = sort(surf.f(i,:));
      ipairs(jaux(i),:) = [aux(1) aux(2)];
      ipairs(jaux(i)+1,:) = [aux(1) aux(3)];
      ipairs(jaux(i)+2,:) = [aux(2) aux(3)];
    endfor
    unique(ipairs,"rows");

    ne = size(ipairs,1);
    ergb = fillrgb(ergb);
    ## pre-allocate edges
    n = rep.nstick;
    rep.nstick += ne;
    [rep ietex] = rep_registertexture(rep,etex);
    empty = struct("name","","x0",[0 0 0],"x1",[0 0 0],"r",erad,...
                   "rgb",ergb,"tex",ietex);
    rep.stick(n+1:rep.nstick) = {empty};
    ## fill edges -> this is the hard part
    for i = n+1:rep.nstick
      rep.stick{i}.x0 = surf.v(ipairs(i-n,1),:);
      rep.stick{i}.x1 = surf.v(ipairs(i-n,2),:);
    endfor
  endif

endfunction

