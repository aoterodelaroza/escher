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

function rep_write_json(obj,file="",nest=0)
% function rep_write_json(obj,file="",nest=0)
%
% rep_write_json - write the representation to a json-style file.
%
% Required input variables:
% rep: representation.
% file: py file. If no file is given, use standard output.
%
% Optional input variables (all have default values):
% {LOG}: verbose level (0=silent,1=verbose).
%

  nn = nest+1;

  ## open file
  if (!isempty(file))
    fid = fopen(file,"w");
  else
    fid = stdout();
  endif

  ## header
  if (strcmp(class(obj),"char"))
    fprintf(fid,"\"%s\"",obj);
  elseif (strcmp(class(obj),"double"))
    if (isscalar(obj))
      if (mod(abs(obj),1) == 0)
        fprintf(fid,"%d",obj);
      elseif (iscomplex(obj))
        fprintf(fid,"(%.16f,%.16f)",real(obj),imag(obj));
      else
        fprintf(fid,"%.16f",obj);
      endif
    elseif (isvector(obj))
      fprintf(fid,"[\n");
      for i = 1:length(obj)
        if (i == length(obj))
          nn = 0;
        endif
        rep_write_json(obj(i),file,nn);
      endfor
      fprintf(fid,"]\n");
    else
      fprintf(fid,"[\n");
      s = size(obj);
      for i = 1:s(1)
        if (i == s(1))
          nn = 0;
        endif
        rep_write_json(reshape(squeeze(obj(i,:)),s(2:end)),file,nn);
      endfor
      fprintf(fid,"]\n");
    endif
  elseif (strcmp(class(obj),"struct"))
    fprintf(fid,"{\n");
    a = fieldnames(obj);
    for i = 1:length(a)
      if (i == length(a))
        nn = 0;
      endif
      fprintf(fid,"\"%s\": ",a{i});
      fi = getfield(obj,a{i});
      rep_write_json(fi,file,nn);
    endfor
    fprintf(fid,"}\n");
  elseif (strcmp(class(obj),"cell"))
    fprintf(fid,"[\n");
    for i = 1:length(obj)
      if (i == length(obj))
        nn = 0;
      endif
      rep_write_json(obj{i},file,nn);
    endfor
    fprintf(fid,"]\n");
  endif

  if (nest > 0)
    fprintf(fid,",\n");
  else
    fprintf(fid,"\n");
  endif

#  fprintf(fid,"  \"name\": \"%s\",\n",rep.name);
#
#  ## balls
#  fprintf(fid,"  \"nball\": %d,\n",rep.nball);
#  fprintf(fid,"  \"ball\": [\n");
#  for i = 1:rep.nball
#    fprintf(fid,"  {\n");
#    fprintf(fid,"    \"name\": \"%s\",\n",rep.ball{i}.name);
#    fprintf(fid,"    \"r\": %.10f,\n",rep.ball{i}.r);
#    fprintf(fid,"    \"rgb\": [%.10f, %.10f, %.10f, %.10f, %.10f],\n",fillrgb(rep.ball{i}.rgb));
#    fprintf(fid,"    \"tex\": %d,\n",rep.ball{i}.tex);
#    fprintf(fid,"    \"x\": [%.10f, %.10f, %.10f]\n",rep.ball{i}.x);
#    if (i == rep.nball)
#      fprintf(fid,"  }\n");
#    else
#      fprintf(fid,"  },\n");
#    endif
#  endfor
#  fprintf(fid,"  ]\n");
#
#  ## sticks
#  fprintf(fid,"  \"nstick\": %d,\n",rep.nstick);
#  fprintf(fid,"  \"stick\": [\n");
#  for i = 1:rep.nstick
#    fprintf(fid,"  {\n");
#    fprintf(fid,"    \"name\": \"%s\",\n",rep.stick{i}.name);
#    fprintf(fid,"    \"x0\": [%.10f, %.10f, %.10f]\n",rep.stick{i}.x0);
#    fprintf(fid,"    \"x1\": [%.10f, %.10f, %.10f]\n",rep.stick{i}.x1);
#    fprintf(fid,"    \"r\": %.10f,\n",rep.stick{i}.r);
#    fprintf(fid,"    \"rgb\": [%.10f, %.10f, %.10f, %.10f, %.10f],\n",fillrgb(rep.stick{i}.rgb));
#    fprintf(fid,"    \"tex\": %d,\n",rep.stick{i}.tex);
#    if (i == rep.nstick)
#      fprintf(fid,"  }\n");
#    else
#      fprintf(fid,"  },\n");
#    endif
#  endfor
#  fprintf(fid,"  ]\n");
#
#  ## vertices
#  fprintf(fid,"  \"nvertex\": %d,\n",rep.nvertex);
#  fprintf(fid,"  \"vertex\": [\n");
#  for i = 1:rep.nvertex
#    fprintf(fid,"  {\n");
#    fprintf(fid,"    \"x\": [%.10f, %.10f, %.10f]\n",rep.triangle{i}.x);
#    fprintf(fid,"    \"tex\": %d,\n",rep.vertex{i}.tex);
#    if (i == rep.nvertex)
#      fprintf(fid,"  }\n");
#    else
#      fprintf(fid,"  },\n");
#    endif
#  endfor
#  fprintf(fid,"  ]\n");
#
#  ## triangles
#  fprintf(fid,"  \"ntriangle\": %d,\n",rep.ntriangle);
#  fprintf(fid,"  \"triangle\": [\n");
#  for i = 1:rep.ntriangle
#    fprintf(fid,"  {\n");
#    fprintf(fid,"    \"idx\": [%d, %d, %d]\n",rep.triangle{i}.idx);
#    fprintf(fid,"    \"tex\": %d,\n",rep.triangle{i}.tex);
#    if (i == rep.ntriangle)
#      fprintf(fid,"  }\n");
#    else
#      fprintf(fid,"  },\n");
#    endif
#  endfor
#  fprintf(fid,"  ]\n");

  if (!isempty(file))
    fclose(fid);
  endif

endfunction
