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

function cr_write_tessel(cr,file="",LOG=0)
% function cr_write_tessel(cr,file="")
%
% cr_write_tessel -- write a basic tessel input file.
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

  fprintf(lu,"set camangle 75 -10 45\n");
  fprintf(lu,"set background background {color rgb <1,1,1>}\n");
  fprintf(lu,"set use_planes .false.\n");
  fprintf(lu,"set ball_texture finish{specular 0.2 roughness 0.1 reflection 0.1}\n");
  fprintf(lu,"set equalscale noscale\n");
  fprintf(lu,"molecule\n");
  fprintf(lu," crystal\n");
  if (isfield(cr,"name"))
    fprintf(lu,"  title %s\n",cr.name);
  endif
  fprintf(lu,"  spg P 1\n");
  fprintf(lu,"  cell %.10f %.10f %.10f %.6f %.6f %.6f\n",cr.a,cr.b*180/pi);
  fprintf(lu,"  crystalbox  -2.30 -2.30 -2.30 2.30 2.30 2.30\n");
  fprintf(lu,"  clippingbox -0.02 -0.02 -0.02 1.02 1.02 1.02\n");
  for i = 1:cr.nat
    fprintf(lu,"  neq   %.10f %.10f %.10f %s\n",cr.x(i,:),cr.attyp{cr.typ(i)});
  endfor
  fprintf(lu," endcrystal\n");
  fprintf(lu," # molmotif allmaincell\n");
  fprintf(lu," # vrml cell.wrl\n");
  fprintf(lu," # povray cell.pov\n");
  fprintf(lu," # vmd cell.vmd\n");
  fprintf(lu,"endmolecule\n");
  fprintf(lu,"# run povray -d +ft +Icell.pov +Ocell.tga +W2000 +H2000 +A\n");
  fprintf(lu,"# run convert cell.tga -bordercolor white -border 1x1 -trim +repage cell.png\n");
  fprintf(lu,"# run rm -f cell.tga\n");
  fprintf(lu,"end\n");

  if (!isempty(file))
    fclose(lu);
  endif

  if (LOG > 0)
    printf("cr_write_tessel: Writing %s\n", file);
    cr_popinfo(cr)
  endif

endfunction
