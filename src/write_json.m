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

function write_json(obj,file="",nest=0)
% function write_json(obj,file="",nest=0)
%
% write_json - write an object to a json-style file.
%
% Required input variables:
% obj: the object to be written. It can be anything: molecules, representations,
%      crystals, etc.
% file: the output json file or standard output (default).
% nest: internal. DO NOT USE.
%

  nn = nest+1;

  ## open file
  if (nest == 0)
    if (!isempty(file))
      fid = fopen(file,"w");
      doclose = 1;
    else
      fid = stdout();
      doclose = 0;
    endif
    file = fid;
  else
    fid = file;
  endif

  ## header
  if (strcmp(class(obj),"char"))
    fprintf(fid,"%*s\"%s\"",nest,"",obj);
  elseif (strcmp(class(obj),"double"))
    if (isscalar(obj))
      if (mod(abs(obj),1) == 0)
        fprintf(fid,"%*s%d",nest,"",obj);
      elseif (iscomplex(obj))
        fprintf(fid,"%*s(%.16f,%.16f)",nest,"",real(obj),imag(obj));
      else
        fprintf(fid,"%*s%.16f",nest,"",obj);
      endif
    elseif (isvector(obj))
      fprintf(fid,"[");
      for i = 1:length(obj)
        write_json(obj(i),file,nn);
        if (i < length(obj))
          fprintf(fid,",");
        endif
      endfor
      fprintf(fid,"]");
    else
      fprintf(fid,"[\n");
      s = size(obj);
      for i = 1:s(1)
        write_json(reshape(squeeze(obj(i,:)),s(2:end)),file,nn);
        if (i < s(1))
          fprintf(fid,",\n");
        endif
      endfor
      fprintf(fid,"]");
    endif
  elseif (strcmp(class(obj),"struct"))
    fprintf(fid,"{\n");
    a = fieldnames(obj);
    for i = 1:length(a)
      fprintf(fid,"%*s\"%s\": ",nest,"",a{i});
      fi = getfield(obj,a{i});
      write_json(fi,file,nn);
      if (i < length(a))
        fprintf(fid,",\n");
      endif
    endfor
    fprintf(fid,"}");
  elseif (strcmp(class(obj),"cell"))
    if (isempty(obj) || isvector(obj))
      fprintf(fid,"[");
      for i = 1:length(obj)
        write_json(obj{i},file,nn);
        if (i < length(obj))
          fprintf(fid,",");
        endif
      endfor
      fprintf(fid,"]");
    else
      fprintf(fid,"[\n");
      s = size(obj);
      for i = 1:s(1)
        write_json(reshape(squeeze(obj(i,:)),s(2:end)),file,nn);
        if (i < s(1))
          fprintf(fid,",\n");
        endif
      endfor
      fprintf(fid,"]");
    endif
  endif

  if (nest == 0)
    fprintf(fid,"\n");
    if (doclose)
      fclose(fid);
    endif
  endif

endfunction
