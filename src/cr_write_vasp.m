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

function cr_write_vasp(cr,file="",LOG=0)
% function cr_write_vasp(cr,file="")
%
% cr_write_vasp -- write a vasp POSCAR file to file or to the stdout if
%                  none is given.
%
% Input:
% cr: crystal structure.
% file: name of the output file. If no file is given (default), then
%       the POSCAR is written to the stdout.
%
% Authors: AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
%          VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
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
  fprintf(lu,"%s (",strtrim(cr.name));
  fprintf(lu,"%s ",cr.attyp{1:cr.ntyp});
  fprintf(lu,")\n");
  fprintf(lu,"%.5f\n",1);
  r = cr.r';
  r = r * bohrtoans;
  fprintf(lu,"%.12f %.12f %.12f\n %.12f %.12f %.12f\n %.12f %.12f %.12f\n",r');
  
  for i = 1:cr.ntyp
    fprintf(lu,"%d ",sum(cr.typ == i));
  endfor
  fprintf(lu,"\n");
  fprintf(lu,"Direct\n");

  for i = 1:cr.ntyp
    for j = 1:cr.nat
      if (cr.typ(j) == i)
        fprintf(lu," %.15f %.15f %.15f \n",cr.x(j,:))
      endif
    endfor
  endfor

  if (!isempty(file))
    fclose(lu);
  endif

  if (LOG > 0)
    printf("cr_write_vasp: Writing %s\n", file);
    cr_popinfo(cr)
  endif

endfunction
