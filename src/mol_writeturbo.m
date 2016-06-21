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

function err = mol_writeturbo(mol, filename="")
%function err = mol_writeturbo(mol, filename="")
%
% mol_writegjf - write a sequence for turbomole's define program
%

  err = 1;
  if (isempty(filename))
     fid = stdout();
  else
     [fid,msg] = fopen(filename,"w");
     if (fid < 0 || ferror(fid)) 
       disp(msg)
       error("mol_readgjf: Could not find -- %s",filename);
     endif
  endif

  fprintf(fid,"\n");
  fprintf(fid,"title\n");
  fprintf(fid,"ai\n");
  for i = 1 : mol.nat;
     fprintf(fid,"%s\n", mol.atname{i});
     fprintf(fid,"%15.9f %15.9f %15.9f A\n", mol.atxyz(1:3,i));
     fprintf(fid,"*\n");
  endfor
  fprintf(fid,"*\n");
  fprintf(fid,"*\n");
  fprintf(fid,"no\n");
  fprintf(fid,"bb all aug-cc-pVTZ\n");
  fprintf(fid,"*\n");
  fprintf(fid,"eht\n");
  fprintf(fid,"y\n");
  fprintf(fid,"0\n");
  fprintf(fid,"y\n");
  fprintf(fid,"dft\n");
  fprintf(fid,"func\n");
  fprintf(fid,"b-lyp\n");
  fprintf(fid,"grid\n");
  fprintf(fid,"m5\n");
  fprintf(fid,"on\n");
  fprintf(fid,"*\n");
  fprintf(fid,"ri\n");
  fprintf(fid,"on\n");
  fprintf(fid,"*\n");
  fprintf(fid,"*\n");

  if (!isempty(filename))
     fclose(fid);
  endif
  err = 0;

endfunction
