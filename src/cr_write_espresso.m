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

function cr_write_espresso(cr,file="")
% function cr_write_espresso(cr,file="")
%
% cr_write_espresso -- write a basic quantum espresso input file.
%
% Input:
% cr: crystal structure.
% file: name of the output file. If no file is given (default), then
%       the file is written to the stdout.
%
% Authors: AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
%          VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
% Created: Nov. 2012

  if (!isstruct(cr) || isempty(cr))
    error("Invalid cr input")
  endif
  if (!isempty(file))
    lu = fopen(file,"w");
  else
    lu = stdout();
  endif

  fprintf(lu,"&control\n");
  fprintf(lu," title='crystal',\n");
  fprintf(lu," prefix='crystal',\n");
  fprintf(lu," pseudo_dir='.',\n");
  fprintf(lu,"/\n");
  fprintf(lu,"&system\n");
  fprintf(lu," ibrav=0,\n");
  fprintf(lu," celldm(1)=1.0,\n");
  fprintf(lu," nat=%d,\n",cr.nat);
  fprintf(lu," ntyp=%d,\n",cr.ntyp);
  fprintf(lu," ecutwfc=60.0,\n");
  fprintf(lu," ecutrho=600.0,\n");
  fprintf(lu," xdm=.true.,\n");
  fprintf(lu," xdm_a1=0.4073,\n");
  fprintf(lu," xdm_a2=2.4150,\n");
  fprintf(lu,"/\n");
  fprintf(lu,"&electrons\n");
  fprintf(lu," conv_thr = 1d-8,\n");
  fprintf(lu,"/\n");
  fprintf(lu,"ATOMIC_SPECIES\n");
  for i = 1:cr.ntyp
    [z,atom] = mol_dbatom(cr.attyp{i});
    fprintf(lu,"%s %14.8f %s.UPF\n",tolower(cr.attyp{i}),atom.mass,tolower(cr.attyp{i}));
  endfor
  fprintf(lu,"\n");
  fprintf(lu,"ATOMIC_POSITIONS crystal\n");
  for i = 1:cr.nat
    fprintf(lu,"%s %.10f %.10f %.10f\n",tolower(cr.attyp{cr.typ(i)}),cr.x(i,:));
  endfor
  fprintf(lu,"\n");
  fprintf(lu,"K_POINTS automatic\n");
  fprintf(lu,"4 4 4 1 1 1\n");
  fprintf(lu,"\n");
  fprintf(lu,"CELL_PARAMETERS cubic\n");
  fprintf(lu," %.10f %.10f %.10f\n",cr.r(1,:));
  fprintf(lu," %.10f %.10f %.10f\n",cr.r(2,:));
  fprintf(lu," %.10f %.10f %.10f\n",cr.r(3,:));

  if (!isempty(file))
    fclose(lu);
  endif

endfunction
