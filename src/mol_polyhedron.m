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

function [rep molc1 molv1]  = mol_polyhedron(molc, molv, addto="", at="", by="", dist=[-1 1.15],...
                                frgb=[0 0 128 115 115], ergb=[0 0 128], ftex="opaque_triangle_default", etex="stick_default",...
                                erad=0.025)
% function [rep molc1 molv1] = mol_polyhedron(molc, molv, rep="", at="", by="", dist=[-1 1.15], \
%                               frgb=[0 0 128 0 225], ergb=[0 0 128], ftex="opaque_triangle_default", etex="stick_default",
%                               erad=0.025)
%
% mol_polyhedron - create polyhedra by giving the center atoms, vertex atoms and a distance criterion
%
% Required input variables:
% molc: molecule containing the centers of the polyhedra.
% molv: molecule containing the vertices of the polyhedra. Usually a superset of the previous.
%       Can be the same as molc.
% at: central atom. Can be:
%     - number: atomic number.
%     - vector: coordinates of an atom (for a single polyhedron)
%     - string: a regexp for the atomic symbol
% by: vertex atom. Can be:
%     - number: atomic number for the vertices
%     - array of numbers: list of atomic numbers for the vertices
%     - string: regexp for the atomic symbols
%     - cell array of strings: regexps for the atomic symbols
%
% Optional input variables (all have default values):
% addto: add the new graphics elements to this representation.
% dist: minimum and maximum atomic distance to place a polyhedron, in the form
%       dist = [dmin dmax]. If the form [-1 f] is used, then place a polyhedron 
%       if the atomic distance is less than f times the sum of covalent radii.
% frgb: rgbft for the face color (from 0 to 255). Make it empty ("" or []) to 
%       deactivate faces.
% ergb: rgbft for the edge color (from 0 to 255). Make it empty ("" or []) to 
%       deactivate edges. 
% ftex: face texture.
% etex: edge texture.
% erad: edge radius.
%
% Output variables:
% rep: the representation containing the polyhedra
% molc1: the molecule containing the polyhedra centers (subset of molc).
% molv1: the molecule containing the polyhedra vertices (subset of molv).
%

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
    if (isfield(molc,"name") && !isempty(molc.name))
      rep.name = molc.name;
    endif
  endif

  ## covalent radii
  nc = molc.nat;
  nv = molv.nat;
  if (dist(1) < 0)
    rcovc = zeros(1,nc);
    for i = 1:nc
      rcovc(i) = mol_rcov(molc.atnumber(i));
    endfor
    rcovv = zeros(1,nv);
    for i = 1:nv
      rcovv(i) = mol_rcov(molv.atnumber(i));
    endfor
  endif

  ## build the polyhedra
  if (!isempty(frgb))
    [rep iftex] = rep_registertexture(rep,ftex);
  endif
  if (!isempty(ergb))
    [rep ietex] = rep_registertexture(rep,etex);
  endif
  ic1 = [];
  iv1 = [];
  for ic = 1:nc
    useit = (ischar(at) && regexp(molc.atname{ic},at));
    useit = useit || (isscalar(at) && molc.atnumber(ic) == at);
    useit = useit || (length(at)==3 && norm(molc.atxyz(:,ic)-at)<1e-10);
    if (useit)
      idx = [];
      for iv = 1:nv
        useit = 0;
        if (iscell(by))
          for k = 1:length(by)
            useit = useit || regexp(molv.atname{iv},by{k});
          endfor
        elseif (ischar(by))
          useit = useit || regexp(molv.atname{iv},by);
        elseif (isscalar(by))
          useit = useit || regexp(molv.atnumber(iv),by);
        else
          for k = 1:length(by)
            useit = useit || (molv.atnumber(iv) == by(k));
          endfor
        endif
        if (!useit)
          continue
        endif
        d = norm(molc.atxyz(:,ic) - molv.atxyz(:,iv));
        if (dist(1) < 0)
          useit = (d <= dist(2) * (rcovv(iv) + rcovc(ic)));
        else
          useit = (d >= dist(1)) && (d <= dist(2));
        endif
        if (useit) 
          idx = [idx iv];
        endif
      endfor

      ## skip polyhedra with less than 4 points
      if (length(idx) < 4)
        continue
      endif

      ic1 = [ic1 ic];
      iv1 = unique([iv1 idx]);

      ## build the polyhedron using the convex hull
      v = molv.atxyz(:,idx)';
      h = convhulln(v);
      if (!isempty(frgb))
        nv0 = rep.nvertex;
        for i = 1:length(idx)
          rep.nvertex += 1;
          rep.vertex{rep.nvertex} = vertex();
          rep.vertex{rep.nvertex}.x = molv.atxyz(:,idx(i))';
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
        icon = zeros(length(idx));
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
    endif
  endfor

  molc1 = mol_getfragment(molc,ic1);
  molv1 = mol_getfragment(molv,iv1);

endfunction
