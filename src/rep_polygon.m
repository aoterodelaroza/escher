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

function rep = rep_polygon(addto="", x0, frgb=[0 0 128 0 128], ergb=[0 0 128], ...
                           ftex="opaque_triangle_default", etex="stick_default", erad=0.025)
% function rep = rep_polygon(addto="", x0, frgb=[0 0 128 0 128], ergb=[0 0 255], ...
%                            ftex="opaque_triangle_default", etex="stick_default", erad=0.025)
%
% rep_polygon - draw a polygon by calculating the barycenter of the vertices in the input
%               and drawing triangles from it to consecutive pairs of vertices.
%               Does *not* check in-plane or convexity.
%
% Required input variables:
% x0: list (nx3) of vertices for the polygon.
%
% Optional input variables (all have default values):
% addto: add the new graphics elements to this representation.
% frgb: rgbft for the face color (from 0 to 255). Make it empty ("" or []) to 
%       deactivate faces.
% ergb: rgbft for the edge color (from 0 to 255). Make it empty ("" or []) to 
%       deactivate edges. 
% ftex: face texture.
% etex: edge texture.
% erad: edge radius.
%

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation_();
    if (isfield(mol,"name") && !isempty(mol.name))
      rep.name = mol.name;
    endif
  endif

  ## center of mass
  xcm = sum(x0,1) / size(x0,1);

  if (!isempty(frgb))
    ## add vertices -> right now not checking a thing
    nv0 = rep.nvertex + 1;
    rep.nvertex += 1; 
    rep.vertex{rep.nvertex} = vertex();
    rep.vertex{rep.nvertex}.x = xcm;
    rep.vertex{rep.nvertex}.rgb = fillrgb(frgb);
    for i = 1:size(x0,1)
      rep.nvertex += 1; 
      rep.vertex{rep.nvertex} = vertex();
      rep.vertex{rep.nvertex}.x = x0(i,:);
      rep.vertex{rep.nvertex}.rgb = fillrgb(frgb);
    endfor

    ## add triangles -> not checking either
    [rep iftex] = rep_registertexture(rep,ftex);
    for i = 1:size(x0,1)
      rep.ntriangle += 1;
      rep.triangle{rep.ntriangle} = triangle();
      rep.triangle{rep.ntriangle}.idx = [nv0, nv0+i, nv0+mod(i,size(x0,1))+1];
      rep.triangle{rep.ntriangle}.rgb = fillrgb(frgb);
      rep.triangle{rep.ntriangle}.tex = iftex;
    endfor
  endif  

#  ## covalent radii
#  nc = molc.nat;
#  nv = molv.nat;
#  if (dist(1) < 0)
#    rcovc = zeros(1,nc);
#    for i = 1:nc
#      [sym atom] = mol_dbsymbol(molc.atnumber(i));
#      rcovc(i) = atom.rcov;
#    endfor
#    rcovv = zeros(1,nv);
#    for i = 1:nv
#      [sym atom] = mol_dbsymbol(molv.atnumber(i));
#      rcovv(i) = atom.rcov;
#    endfor
#  endif
#
#  for ic = 1:nc
#    useit = (ischar(at) && regexp(molc.atname{ic},at));
#    useit = useit || (isscalar(at) && molc.atnumber(ic) == at);
#    useit = useit || (length(at)==3 && norm(molc.atxyz(:,ic)-at)<1e-10);
#    if (useit)
#      idx = [];
#      for iv = 1:nv
#        useit = 0;
#        if (iscell(by))
#          for k = 1:length(by)
#            useit = useit || regexp(molv.atname{iv},by{k});
#          endfor
#        elseif (ischar(by))
#          useit = useit || regexp(molv.atname{iv},by);
#        elseif (isscalar(by))
#          useit = useit || regexp(molv.atnumber(iv),by);
#        else
#          for k = 1:length(by)
#            useit = useit || (molv.atnumber(iv) == by(k));
#          endfor
#        endif
#        if (!useit)
#          continue
#        endif
#        d = norm(molc.atxyz(:,ic) - molv.atxyz(:,iv));
#        if (dist(1) < 0)
#          useit = (d <= dist(2) * (rcovv(iv) + rcovc(ic)));
#        else
#          useit = (d >= dist(1)) && (d <= dist(2));
#        endif
#        if (useit) 
#          idx = [idx iv];
#        endif
#      endfor
#
#      ## build the polyhedron using the convex hull
#      v = molv.atxyz(:,idx)';
#      h = convhulln(v);
#      if (!isempty(frgb))
#        nv0 = rep.nvertex;
#        for i = 1:length(idx)
#          rep.nvertex += 1;
#          rep.vertex{rep.nvertex}.x = molv.atxyz(:,idx(i))';
#          rep.vertex{rep.nvertex}.rgb = fillrgb(frgb);
#        endfor
#        nv0 = rep.ntriangle;
#        [rep iftex] = rep_registertexture(rep,ftex);
#        for i = 1:size(h,1)
#          rep.ntriangle += 1;
#          rep.triangle{rep.ntriangle}.idx = nv0 + h(i,:);
#          rep.triangle{rep.ntriangle}.rgb = fillrgb(frgb);
#          rep.triangle{rep.ntriangle}.tex = iftex;
#        endfor
#      endif
#      if (!isempty(ergb))
#        icon = zeros(length(idx));
#        kk = [1 2; 1 3; 2 3];
#        [rep ietex] = rep_registertexture(rep,etex);
#        for i = 1:size(h,1)
#          for j = 1:3
#            rep.nstick = rep.nstick + 1;
#            rep.stick{rep.nstick}.name = "";
#            rep.stick{rep.nstick}.x0 = v(h(i,kk(j,1)),:);
#            rep.stick{rep.nstick}.x1 = v(h(i,kk(j,2)),:);
#            rep.stick{rep.nstick}.r = erad;
#            rep.stick{rep.nstick}.rgb = ergb;
#            rep.stick{rep.nstick}.tex = ietex;
#          endfor
#        endfor
#      endif
#    endif
#  endfor

endfunction
