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

function [cr] = cr_read_vasp(file="CONTCAR",potcar="POTCAR",LOG=0)
% function [cr] = cr_readvasp(file="CONTCAR",potcar="POTCAR",LOG=0)
%
% cr_readvasp -- read a VASP crystal structure (POSCAR or CONTCAR)
%                and the atomic types from the POTCAR (or from a cell 
%                array provided by the user).
%
% Input:
% file: the name of the POSCAR or CONTCAR file.
% potcar: the name of the POTCAR file or, alternatively, a cell array
%         containing the atomic symbols (e.g. {"Si","O"}).
%
% Output: 
% cr: crystal description. If no POTCAR is provided, the atom types
%     will be missing from the structure.
%
% Authors: AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
%          VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
% Created: Nov. 2012

  bohrtoang = 0.52917720859;

  cr = struct();
  if (!exist(file,"file"))
    error(sprintf("Could not find file: %s\n",file));
  endif

  ## open file and read title
  fid = fopen(file,"r");
  cr.name = fgetl(fid);
  
  ## cell
  fac = fscanf(fid,"%f",1);
  r = fscanf(fid,"%f",9);
  r = reshape(r,3,3)';
  if (fac < 0) 
    fac = (abs(fac / det(r)))^(1/3);
  endif
  r = r * fac / bohrtoang;

  cr.r = r;
  cr.g = cr.r * cr.r';
  cr.a = sqrt(diag(cr.g))';
  cr.b(1) = acos(cr.g(2,3) / (cr.a(2)*cr.a(3)));
  cr.b(2) = acos(cr.g(1,3) / (cr.a(1)*cr.a(3)));
  cr.b(3) = acos(cr.g(1,2) / (cr.a(1)*cr.a(2)));
  cr.omega = det(cr.r);

  ## read atom types: 2 types of POSCAR, ignore the atomic symbols
  nat = fscanf(fid,"%d"); 
  if (isempty(nat))
    line = fgetl(fid);
    nat = fscanf(fid,"%d");
  endif  
  nat = nat';
  tat = sum(nat);

  cr.ntyp = length(nat);
  typcount = nat;
  cr.nat = tat;
  cr.typ = zeros(1,tat);
  nn = 0;
  for i = 1:cr.ntyp
    for j = 1:typcount(i)
      nn += 1;
      cr.typ(nn) = i;
    endfor
  endfor

  ## ignore direct/selective dynamics line
  dum = fgetl(fid);
  if (strcmpi(substr(dum,1,1),"s"))
    dum = fgetl(fid);
  endif
  iscart = 0;
  if (strcmpi(substr(dum,1,1),"d"))
    iscart = 0;
  elseif (strcmpi(substr(dum,1,1),"c") || strcmpi(substr(dum,1,1),"k"))
    iscart = 1;
    ri = inv(r * bohrtoang);
  else
    error("Unknown cartesian/direct option");
  endif

  ## read crystallographic coordinates
  cr.x = zeros(tat,3);
  for i = 1:tat
    line = fgetl(fid);
    cr.x(i,1:3) = sscanf(line,"%f %f %f");
    if (iscart) 
      cr.x(i,1:3) = cr.x(i,1:3) * ri;
    endif
  endfor
  fclose(fid);

  ## read POTCAR
  dum = cell(10,1);
  cr.ztyp = zeros(1,cr.ntyp);
  cr.zvaltyp = zeros(1,cr.ntyp);
  cr.attyp = cell(1,cr.ntyp);
  for i = 1:length(cr.attyp)
    cr.attyp{i} = "Bq";
  endfor
  if (ischar(potcar) && exist(potcar,"file"))
    fid = fopen(potcar,"r");

    ## read it
    nn = 0;
    do 
      line = fgetl(fid);
      word = sscanf(line,"%s","C");
      if (strcmp(word,"TITEL"))
        nn += 1;
        [dum{1:3} word] = sscanf(line,"%s %s %s %s","C");
        cr.attyp{nn} = word;
      elseif (strcmp(word,"POMASS"))
        [dum{1:5} cr.zvaltyp(nn)] = sscanf(line,"%s %s %s %s %s %f","C");
      endif
    until (feof(fid))
    fclose(fid);
  elseif (iscell(potcar))
    if (length(potcar) != cr.ntyp) 
      error("Wrong number of atom types in cell array");
    endif
    for i = 1:length(potcar)
      cr.attyp{i} = potcar{i};
    endfor
  else
    error(sprintf("Not a cell array and could not find file: %s\n",potcar));
  endif

  ## convert to Z
  for i = 1:cr.ntyp
    cr.ztyp(i) = mol_dbatom(cr.attyp{i},0);
  endfor

  if (LOG > 0)
    printf("cr_read_vasp: Reading %s\n", file);
    cr_popinfo(cr)
  endif

endfunction

