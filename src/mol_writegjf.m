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

function err = mol_writegjf(mol, filename="none")
%function err = mol_writegjf(mol, filename="none")
%
% mol_writegjf - write a gaussian09 input file.
%

  err = 1;
  if (strcmpi(filename,"none"))
     fid = stdout();
  else
     [fid,msg] = fopen(filename,"w");
     if (fid < 0 || ferror(fid)) 
       disp(msg)
       error("mol_writegjf: Could not open -- %s",filename);
     endif
  endif

  fprintf(fid, "%%mem=2GB\n");
  fprintf(fid, "%%nproc=4\n");
  fprintf(fid, "#p b3lyp sto-3g\n");
  fprintf(fid, "\n");
  fprintf(fid, "title\n");
  fprintf(fid, "\n");
  fprintf(fid, "0 1\n");
  for i = 1 : mol.nat;
     fprintf(fid,"   %-2s  ", mol.atname{i});
     fprintf(fid,"  %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
     fprintf(fid,"\n");
  endfor
  fprintf(fid, "\n");

  if (!strcmpi(filename,"none"))
     fclose(fid);
  endif
  err = 0;

endfunction
