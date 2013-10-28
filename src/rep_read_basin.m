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

function rep = rep_read_basin(file, addto="", frgb=[255 128 128 0 128], ergb=[0 0 128], \
                              ftex="opaque_triangle_default", etex="stick_default", erad=0.005)
% function rep = rep_read_basin(file, addto="", frgb=[255 128 128 0 128], ergb=[0 0 128], \
%                               ftex="opaque_triangle_default", etex="stick_default", erad=0.005)
%
% rep_read_basin - read a basin file and create a surface representation.
%
% Required input variables:
% file: the basin file
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

  bohrtoans = 0.52917720859;

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
    rep.name = "basin";
  endif

  ## open the file
  if (!exist(file,"file"))
    error(sprintf("Could not find file: %s\n",file));
  endif
  fid = fopen(file,"r");
  do 
    line = fgetl(fid);
    [first,count] = sscanf(line," %1c",1);
  until (first != "#" && first != " ")
  [nv,nf,ne] = sscanf(line,"%d %d %d","C");
  line = fgetl(fid);
  line = fgetl(fid);

  ## read the vertices
  xv = zeros(nv,3);
  for i = 1:nv
    line = fgetl(fid);
    xv(i,:) = sscanf(line,"%f",3);
  endfor
  xv = xv * bohrtoans;

  ## fill the colors
  if (!isempty(ergb))
    ergb = fillrgb(ergb);
  endif
  if (!isempty(frgb))
    frgb = fillrgb(frgb);
  endif

  ## fill the vertices
  if (!isempty(frgb))
    nv0 = rep.nvertex;
    for i = 1:nv
      rep.nvertex += 1;
      rep.vertex{rep.nvertex} = vertex();
      rep.vertex{rep.nvertex}.x = xv(i,:);
      rep.vertex{rep.nvertex}.rgb = frgb;
    endfor
  endif
  
  ## read the faces
  for i = 1:nf
    line = fgetl(fid);
    idx = sscanf(line,"%d",5);
    idx = idx(2:end) + 1;
    ## edges
    if (!isempty(ergb))
      if (length(idx) == 3)
        kk = [1 2; 1 3; 2 3];
      elseif (length(idx) == 4)
        kk = [1 2; 1 3; 2 3; 1 4; 3 4];
      else
        error("don't know how to handle these polygons");
      endif
      for j = 1:size(kk,1)
        rep.nstick++;
        rep.stick{rep.nstick} = stick();
        rep.stick{rep.nstick}.name = "";
        rep.stick{rep.nstick}.x0 = xv(idx(kk(j,1)),:);
        rep.stick{rep.nstick}.x1 = xv(idx(kk(j,2)),:);
        rep.stick{rep.nstick}.r = erad;
        rep.stick{rep.nstick}.rgb = ergb;
        rep.stick{rep.nstick}.tex = etex;
      endfor
    endif
    ## edges
    if (!isempty(frgb))
      if (length(idx) == 3)
        rep.ntriangle += 1;
        rep.triangle{rep.ntriangle} = triangle();
        rep.triangle{rep.ntriangle}.idx = nv0 + [idx(1) idx(2) idx(3)];
        rep.triangle{rep.ntriangle}.rgb = frgb;
        rep.triangle{rep.ntriangle}.tex = ftex;
      elseif (length(idx) == 4)
        rep.ntriangle += 1;
        rep.triangle{rep.ntriangle} = triangle();
        rep.triangle{rep.ntriangle}.idx = nv0 + [idx(1) idx(2) idx(3)];
        rep.triangle{rep.ntriangle}.rgb = frgb;
        rep.triangle{rep.ntriangle}.tex = ftex;
        rep.ntriangle += 1;
        rep.triangle{rep.ntriangle} = triangle();
        rep.triangle{rep.ntriangle}.idx = nv0 + [idx(1) idx(3) idx(4)];
        rep.triangle{rep.ntriangle}.rgb = frgb;
        rep.triangle{rep.ntriangle}.tex = ftex;
      else
        error("don't know how to handle these polygons");
      endif
    endif
  endfor

  fclose(fid);

endfunction
