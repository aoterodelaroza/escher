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

function err = mol_writeturbo(mol, filename="", bsic="")
%function err = mol_writeturbo(mol, filename="", bsic="")
%
% mol_writegjf - write a sequence for turbomole's define program.
%   If filename is not given, use stdout. If bsic is non-empty,
%   use the tricks for BSIC calculations. Pass the name of the 
%   basis set file in the bsic keyword (e.g. minis).
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

  ats = cell(mol.nat,1);
  atz = zeros(mol.nat,1);
  for i = 1:mol.nat
    atz(i) = mol_dbatom(mol.atname{i});
    ats{i} = lower(mol_dbsymbol(atz(i)));
  endfor
  uatz = unique(atz);
  uats = cell(length(uatz),1);
  for i = 1:length(uatz)
    uats{i} = lower(mol_dbsymbol(uatz(i)));
  endfor

  fprintf(fid,"\n");
  fprintf(fid,"title\n");
  fprintf(fid,"ai\n");
  for i = 1:mol.nat;
     fprintf(fid,"%s\n", ats{i});
     fprintf(fid,"%15.9f %15.9f %15.9f A\n", mol.atxyz(1:3,i));
     fprintf(fid,"*\n");
  endfor
  fprintf(fid,"*\n");
  fprintf(fid,"*\n");
  fprintf(fid,"no\n");
  if (!isempty(bsic))
    for i = 1:length(uats)
      fprintf(fid,"b \"%s\" %s\n",uats{i},bsic);
      fprintf(fid,"file %s.tbs\n",bsic);
    endfor
  else
    fprintf(fid,"bb all aug-cc-pVTZ\n");
  endif
  fprintf(fid,"*\n");
  fprintf(fid,"eht\n");
  fprintf(fid,"y\n");
  fprintf(fid,"0\n");
  fprintf(fid,"y\n");
  if (!isempty(bsic)) 
    fprintf(fid,"&\n");
    fprintf(fid,"&\n");
    fprintf(fid,"y\n");
    for i = 1:length(uats)
      fprintf(fid,"c \"%s\" %d\n",uats{i},uatz(i)+2);
    endfor
    for i = 1:length(uats)
      fprintf(fid,"ecp \"%s\" %s-bsic\n",uats{i},bsic);
      fprintf(fid,"file %s.tbsic\n",bsic);
    endfor
    fprintf(fid,"*\n");
    fprintf(fid,"*\n");
    fprintf(fid,"\n");
  endif
  fprintf(fid,"dft\n");
  fprintf(fid,"func\n");
  fprintf(fid,"b-lyp\n");
  fprintf(fid,"grid\n");
  fprintf(fid,"m5\n");
  fprintf(fid,"on\n");
  fprintf(fid,"*\n");
  fprintf(fid,"ri\n");
  fprintf(fid,"off\n");
  fprintf(fid,"*\n");
  fprintf(fid,"*\n");

  if (!isempty(filename))
     fclose(fid);
  endif
  err = 0;

endfunction
