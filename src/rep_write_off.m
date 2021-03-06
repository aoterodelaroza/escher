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

function rep_write_off(rep,file="",ndiv=1,nbase=16,LOG=0)
% function rep_write_off(rep,file="",ndiv=1,nbase=16,LOG=0)
%
% rep_write_off - write a representation to a geomview (off) file.
%
% Required input variables:
% rep: representation.
% file: off file. If no file is given, use standard output.
% ndiv: spheres are represented as icosahedra with ndiv refinmenents.
% nbase: sticks are represented as tubes whose bases are polygons with
%        nbase sides.
%
% Optional input variables (all have default values):
% {LOG}: verbose level (0=silent,1=verbose).
%

  bohrtoans = 0.52917720859;

  ## create the sphere and the cylinder model
  nvsph = nfsph = nvstk = nfstk = 0;
  if (rep.nball > 0)
    [vsph ifsph] = spheremodel(ndiv);
    nvsph = size(vsph,1);
    nfsph = size(ifsph,1);
  endif
  if (rep.nstick > 0)
    [vstk ifstk] = cylindermodel(nbase);
    nvstk = size(vstk,1);
    nfstk = size(ifstk,1);
  endif

  ## open file
  if (!isstruct(rep) || isempty(rep))
    error("Invalid rep input")
  endif
  if (!isempty(file))
    fid = fopen(file,"w");
  else
    fid = stdout();
  endif

  ## count the total number of vertices and faces
  nv = 0; nf = 0;
  nv = nv + nvsph * rep.nball;
  nf = nf + nfsph * rep.nball;

  nn = 0;
  for i = 1:rep.nstick
    rr = rep.stick{i}.x1 - rep.stick{i}.x0;
    if (any(abs(rr) > eps))
      nn++;
    endif
  endfor
  nv = nv + nvstk * nn;
  nf = nf + nfstk * nn;
  nv = nv + rep.nvertex;
  nf = nf + rep.ntriangle;
  
  ## write the polyhedra
  fprintf(fid,"# OFF created by molware\n");
  if (isfield(rep,"name"))
    fprintf(fid,"# title: %s\n",rep.name);
  endif

  ## total number of vertices and faces
  fprintf(fid,"OFF\n");
  fprintf(fid,"%d %d %d\n",nv,nf,0);

  ## ball and stick vertices
  nwri = 0;

  npad_vert = nwri;
  for i = 1:rep.nvertex
    fprintf(fid,"%.9f %.9f %.9f\n",rep.vertex{i}.x);
  endfor
  nwri = nwri + rep.nvertex;

  npad_bal = nwri;
  for i = 1:rep.nball
    for j = 1:nvsph
      fprintf(fid,"%.9f %.9f %.9f\n",rep.ball{i}.x+rep.ball{i}.r*vsph(j,:));
    endfor
  endfor
  nwri = nwri + rep.nball * nvsph;

  npad_stk = nwri;
  for i = 1:rep.nstick
    ## local coordinates of the tube
    rr = rep.stick{i}.x1 - rep.stick{i}.x0;
    if (abs(rr(1)) > eps)
      v1 = [-rr(2)/rr(1), 1, 0];
    elseif (abs(rr(2)) > eps);
      v1 = [0, -rr(3)/rr(2), 1];
    elseif (abs(rr(3)) > eps);
      v1 = [1, 0, -rr(1)/rr(3)];
    else
      ## zero-length stick, skip
      continue
    endif
    v2 = cross(rr,v1);
    v1 = v1 / norm(v1) * rep.stick{i}.r; 
    v2 = v2 / norm(v2) * rep.stick{i}.r;
    
    ## cylinder scaled using the local coordinates
    for j = 1:nvstk
      fprintf(fid,"%.9f %.9f %.9f\n",rep.stick{i}.x0+vstk(j,1)*v1+vstk(j,2)*v2+vstk(j,3)*rr);
    endfor
  endfor
  nwri = nwri + rep.nstick * nvstk;

  ## faces
  for i = 1:rep.ntriangle
    fprintf(fid,"%d %d %d %d\n",3,npad_vert+rep.triangle{i}.idx-1);
  endfor

  for i = 1:rep.nball
    nadd = (i-1)*nvsph;
    for j = 1:nfsph
      fprintf(fid,"%d %d %d %d\n",3,npad_bal+ifsph(j,:)+nadd-1);
    endfor
  endfor

  for i = 1:rep.nstick
    nadd = (i-1)*nvstk;
    for j = 1:nfstk
      fprintf(fid,"%d %d %d %d\n",3,npad_stk+ifstk(j,:)+nadd-1);
    endfor
  endfor

  if (!isempty(file))
    fclose(fid);
  endif
  if (LOG > 0)
    printf("rep_write_off: Writing %s\n", file);
  endif

endfunction
