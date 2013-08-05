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

function rep = mol_polyhedron(molc, molv, addto="", at="", by="", dist=[-1 1.15], \
                              frgb=[0 0 128 0 225], ergb=[0 0 128], ftex="opaque_triangle_default", etex="stick_default",\
                              erad=0.025)
% function rep = mol_polyhedron(molc, molv, rep="", at="", by="", dist=[-1 1.15], \
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

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
    if (isfield(mol,"name") && !isempty(mol.name))
      rep.name = mol.name;
    endif
  endif

  ## initialize
  if (!isfield(rep,"nstick"))
    rep.nstick = 0;
  endif
  if (rep.nstick == 0)
    rep.stick = cell();
  endif
  if (!isfield(rep,"ntriangle"))
    rep.ntriangle = 0;
  endif
  if (rep.ntriangle == 0)
    rep.triangle = cell();
  endif
  if (!isfield(rep,"nvertex"))
    rep.nvertex = 0;
  endif
  if (rep.nvertex == 0)
    rep.vertex = cell();
  endif

  ## covalent radii
  nc = length(molc.atnumber);
  nv = length(molv.atnumber);
  if (dist(1) < 0)
    rcovc = zeros(1,nc);
    for i = 1:nc
      [sym atom] = mol_dbsymbol(molc.atnumber(i));
      rcovc(i) = atom.rcov;
    endfor
    rcovv = zeros(1,nv);
    for i = 1:nv
      [sym atom] = mol_dbsymbol(molv.atnumber(i));
      rcovv(i) = atom.rcov;
    endfor
  endif

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

      ## build the polyhedron using the convex hull
      v = molv.atxyz(:,idx)';
      h = convhulln(v);
      if (!isempty(frgb))
        nv0 = rep.nvertex;
        for i = 1:length(idx)
          rep.nvertex += 1;
          rep.vertex{rep.nvertex}.x = molv.atxyz(:,idx(i))';
          rep.vertex{rep.nvertex}.rgb = fillrgb(frgb);
        endfor
        nv0 = rep.ntriangle;
        for i = 1:size(h,1)
          rep.ntriangle += 1;
          rep.triangle{rep.ntriangle}.idx = nv0 + h(i,:);
          rep.triangle{rep.ntriangle}.rgb = fillrgb(frgb);
          rep.triangle{rep.ntriangle}.tex = ftex;
        endfor
      endif
      if (!isempty(ergb))
        icon = zeros(length(idx));
        kk = [1 2; 1 3; 2 3];
        for i = 1:size(h,1)
          for j = 1:3
            rep.nstick = rep.nstick + 1;
            rep.stick{rep.nstick}.name = "";
            rep.stick{rep.nstick}.x0 = v(h(i,kk(j,1)),:);
            rep.stick{rep.nstick}.x1 = v(h(i,kk(j,2)),:);
            rep.stick{rep.nstick}.r = erad;
            rep.stick{rep.nstick}.rgb = ergb;
            rep.stick{rep.nstick}.tex = etex;
          endfor
        endfor
      endif
    endif
  endfor

endfunction
