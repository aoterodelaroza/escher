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

function rep_write_py(rep,file="")
% function rep_write_py(rep,file="")
%
% rep_write_py - xxxx
%
% Required input variables:
% rep: representation.
% file: py file. If no file is given, use standard output.
%
% Optional input variables (all have default values):
% {LOG}: verbose level (0=silent,1=verbose).
%

  
  ## open file
  if (!isstruct(rep) || isempty(rep))
    error("Invalid rep input")
  endif
  if (!isempty(file))
    fid = fopen(file,"w");
  else
    fid = stdout();
  endif

  fprintf(fid,"{\n");
  fprintf(fid,"  \"name\": \"%s\",\n",rep.name);
  fprintf(fid,"  \"nball\": %d,\n",rep.nball);
  fprintf(fid,"  \"ball\": [\n");
  for i = 1:rep.nball
    fprintf(fid,"  {\n");
    fprintf(fid,"    \"name\": \"%s\",\n",rep.ball{i}.name);
    fprintf(fid,"    \"r\": %.10f,\n",rep.ball{i}.r);
    fprintf(fid,"    \"rgb\": [%.10f, %.10f, %.10f, %.10f, %.10f],\n",fillrgb(rep.ball{i}.rgb));
    fprintf(fid,"    \"tex\": %d,\n",rep.ball{i}.tex);
    fprintf(fid,"    \"x\": [%.10f, %.10f, %.10f]\n",rep.ball{i}.x);
    if (i == rep.nball)
      fprintf(fid,"  }\n");
    else
      fprintf(fid,"  },\n");
    endif
  endfor
  fprintf(fid,"  ]\n");
  fprintf(fid,"}\n");

  if (!isempty(file))
    fclose(fid);
  endif

endfunction
