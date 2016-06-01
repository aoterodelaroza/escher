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

function rep_write_obj(rep,file="",filemtl="",ndiv=1,nbase=16,LOG=0)
% function rep_write_obj(rep,file="",filemtl="",ndiv=1,nbase=16,LOG=0)
%
% rep_write_obj - write a representation to a wavefront OBJ file. 
%
% Required input variables:
% rep: representation.
% file: obj file. If no file is given, use standard output.
% filemtl: mtl file. If no obj file is given, it is not written.
%          If obj is gievn but no mtl file is present, then use the same
%          root as the obj file and .mtl extension.
% ndiv: spheres are represented as icosahedra with ndiv refinmenents.
% nbase: sticks are represented as tubes whose bases are polygons with
%        nbase sides.
%
% Optional input variables (all have default values):
% {LOG}: verbose level (0=silent,1=verbose).
%

  bohrtoans = 0.52917720859;

  ## textures for colors: run over all objects and collect colors in as few textures 
  ## as possible
  nrgb = 0;
  rgb = zeros(1,5);
  tex = cell();

  for i = 1:rep.nball
    rep.ball{i}.rgb = fillrgb(rep.ball{i}.rgb);
    found = 0;
    for j = 1:nrgb
      if (all(abs(rep.ball{i}.rgb - rgb(j,:)) < 1e-5) && strcmpi(tex{j},rep.ball{i}.tex))
        found = 1;
        rep.ball{i}.irgb = j;
      endif
    endfor
    if (!found)
      nrgb = nrgb + 1;
      rgb(nrgb,1:5) = rep.ball{i}.rgb;
      tex{nrgb} = rep.ball{i}.tex;
      rep.ball{i}.irgb = nrgb;
    endif
  endfor

  for i = 1:rep.nstick
    rep.stick{i}.rgb = fillrgb(rep.stick{i}.rgb);
    found = 0;
    for j = 1:nrgb
      if (all(abs(rep.stick{i}.rgb - rgb(j,:)) < 1e-5) && strcmpi(tex{j},rep.stick{i}.tex))
        found = 1;
        rep.stick{i}.irgb = j;
      endif
    endfor
    if (!found)
      nrgb = nrgb + 1;
      rgb(nrgb,1:5) = rep.stick{i}.rgb;
      tex{nrgb} = rep.stick{i}.tex;
      rep.stick{i}.irgb = nrgb;
    endif
  endfor

  ## write the mtl file
  if (!isempty(file))
    if (isempty(filemtl))
      filemtl = strrep(file,".obj",".mtl");
    endif
    fid = fopen(filemtl,"w");
    for i = 1:nrgb
      atex = texture("obj",tex{i});
      fprintf(fid,"newmtl mat%d\n",i);
      fprintf(fid,"Ns  %.5f\n",atex.Ns);
      fprintf(fid,"Ka  %.5f %.5f %.5f\n",atex.Ka/255);
      fprintf(fid,"Kd  %.5f %.5f %.5f\n",rgb(i,1:3)/255);
      fprintf(fid,"Ks  %.5f %.5f %.5f\n",atex.Ks/255);
      fprintf(fid,"Ni  %.5f\n",atex.Ni);
      fprintf(fid,"d  %.5f\n",1-rgb(i,5)/255);
      fprintf(fid,"Tr  %.5f\n",1-rgb(i,5)/255);
      fprintf(fid,"illum %d\n",atex.illum);
      fprintf(fid,"\n\n");
    endfor
    fclose(fid);
  endif

  ## Create the sphere and the cylinder
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

  ## write the polyhedra
  fprintf(fid,"# OBJ created by molware\n");
  if (isfield(rep,"name"))
    fprintf(fid,"# title: %s\n",rep.name);
  endif
  if (!isempty(file))
    fprintf(fid,"mtllib %s\n",filemtl);
  endif

  ## balls
  nadd = 0;
  for i = 1:rep.nball
    if (isfield(rep.ball{i},"name") && !isempty(rep.ball{i}.name))
      fprintf(fid,"o ball_%s\n",rep.ball{i}.name);
    else
      fprintf(fid,"o ball_%d\n",i);
    endif
    fprintf(fid,"s on\n");
    if (!isempty(file))
      fprintf(fid,"usemtl mat%d\n",rep.ball{i}.irgb);
    endif
    for j = 1:nvsph
      fprintf(fid,"v %.9f %.9f %.9f\n",rep.ball{i}.x+rep.ball{i}.r*vsph(j,:));
    endfor
    for j = 1:nfsph
      fprintf(fid,"f %d %d %d\n",nadd+ifsph(j,:));
    endfor
    nadd = nadd + nvsph;
  endfor

  ## sticks
  for i = 1:rep.nstick
    ## skip degenerate cylinders
    if (norm(rep.stick{i}.x1 - rep.stick{i}.x0) < eps)
      continue
    endif
    ## header
    if (isfield(rep.stick{i},"name") && !isempty(rep.stick{i}.name))
      fprintf(fid,"o stick_%s\n",rep.stick{i}.name);
    else
      fprintf(fid,"o stick_%d\n",i);
    endif
    fprintf(fid,"s on\n");
    if (!isempty(file))
      fprintf(fid,"usemtl mat%d\n",rep.stick{i}.irgb);
    endif

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
      fprintf(fid,"v %.9f %.9f %.9f\n",rep.stick{i}.x0+vstk(j,1)*v1+vstk(j,2)*v2+vstk(j,3)*rr);
    endfor
    for j = 1:nfstk
      fprintf(fid,"f %d %d %d\n",ifstk(j,:)+nadd);
    endfor
    nadd = nadd + nvstk;
  endfor

  if (!isempty(file))
    fclose(fid);
  endif
  if (LOG > 0)
    printf("rep_write_obj: Writing %s\n", file);
  endif

endfunction
