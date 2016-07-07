% Copyright (C) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function cr_write_critic2(cr,file="",LOG=0)
% function cr_write_critic2(cr,file="")
%
% cr_write_critic2 -- write a basic critic2 input file.
%
% Input:
% cr: crystal structure.
% file: name of the output file. If no file is given (default), then
%       the file is written to the stdout.
%
% Authors: AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
%          VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
% Created: Nov. 2012

  bohrtoans = 0.52917720859;

  if (!isstruct(cr) || isempty(cr))
    error("Invalid cr input")
  endif
  if (!isempty(file))
    lu = fopen(file,"w");
  else
    lu = stdout();
  endif

  fprintf(lu,"crystal\n");
  fprintf(lu,"  cell %.10f %.10f %.10f %.6f %.6f %.6f\n",cr.a,cr.b*180/pi);
  for i = 1:cr.nat
    fprintf(lu,"  neq   %.10f %.10f %.10f %s\n",cr.x(i,:),cr.attyp{cr.typ(i)});
  endfor
  fprintf(lu,"endcrystal\n");
  fprintf(lu,"#auto\n");
  fprintf(lu,"#qtree 6\n");
  fprintf(lu,"#nciplot\n");
  fprintf(lu,"#endnciplot\n");
  fprintf(lu,"end\n");

  if (!isempty(file))
    fclose(lu);
  endif

  if (LOG > 0)
    printf("cr_write_critic2: Writing %s\n", file);
    cr_popinfo(cr)
  endif

endfunction
