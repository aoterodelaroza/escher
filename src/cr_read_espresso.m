% Copyright (c) 2012 Alberto Otero-de-la-Roza and Victor Lua~na
% Adapted from an AOR (2011) routine.
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

function cr = cr_read_espresso(file, LOG=0)
% function cr = cr_read_espresso(file,LOG=0)
%
% cr_read_espresso - read in a crystal description from an optimized
% quantum espresso (qe, pwscf) structure.
%
% Required input variables:
% {file}: name of the pwscf output file.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Ouput variables:
% {cr}: struct with the crystals description.
% The main content of this struct includes:
% .name - name of the crystal.
% .nat - number of atoms in the main cell.
% .ntyp - number of atom types.
% .ztyp[1:ntyp] - atomic number for each type.
% .zvaltyp[1:ntyp] - ?
% .attyp{1:ntyp} - atomic symbols for each type.
% .typ[1:nat] - type index for each atom.
% .x[1:nat,1:3] - crystallographic coordinates.
% .r[1:3,1:3] - convert crystallographic to cartesian coordinates.
% .g[1:3,1:3] - crystallographic G matrix.
% .a[1:3] - cell parameters (a,b,c) in bohr (???)
% .b[1:3] - cell angles (alpha,beta,gamm) in radians.
% .omega - cell volume.
%

  bohrtoans = 0.52917720859;

  cr = struct();

  ## open the file
  if (!exist(file,"file"))
    error(sprintf("Could not find file: %s\n",file));
  endif
  fid = fopen(file,"r");

  ## parse the output file
  idocart = -1;
  line = fgetl(fid);
  nocrystalaxes = 0;
  noatpos = 0;
  do 
    ## at the beginning of the run
    if(strfind(line,"Title:"))
      cr.name = fgetl(fid);
    elseif(strfind(line,"bravais-lattice index"))
      ibrav = sscanf(substr(line,index(line,"=")+1),"%d");
    elseif(strfind(line,"lattice parameter (alat)"))
      alat = sscanf(substr(line,index(line,"=")+1),"%d");
    elseif(strfind(line,"number of atoms/cell"))
      cr.nat = sscanf(substr(line,index(line,"=")+1),"%d");
    elseif(strfind(line,"number of atomic types"))
      cr.ntyp = sscanf(substr(line,index(line,"=")+1),"%d");
    elseif(strfind(line,"crystal axes:") && !nocrystalaxes)
      r = zeros(3,3);
      for i = 1:3
        line = fgetl(fid);
        line = substr(line,index(line,"=")+1);
        line = substr(line,index(line,"(")+1);
        r(i,:) = sscanf(line,"%f %f %f");
      endfor
      r = r * alat;
      rinv = inv(r);
    endif

    ## atomic species
    if(regexp(line,"^ *atomic species *valence *mass *pseudopotential *$"))
      cr.zvaltyp = cr.ztyp = zeros(1,cr.ntyp);
      cr.attyp = cell(1,cr.ntyp);
      cr.qtyp = zeros(1,cr.ntyp);
      for i = 1:cr.ntyp
        line = fgetl(fid);
        [dum1 cr.zvaltyp(i) dum2 cr.attyp{i}] = sscanf(line,"%s %f %f %s","C");
        cr.ztyp(i) = mol_dbatom(cr.attyp{i},0);
        cr.attyp{i} = dum1;
      endfor
    endif

    ## coordinates in the input section
    if(regexp(line,"^ *Cartesian axes *$") && !noatpos);
      line = fgetl(fid); line = fgetl(fid);
      cr.typ = zeros(1,cr.nat);
      cr.x = zeros(cr.nat,3);
      for i = 1:cr.nat
        line = fgetl(fid);
        [dum atsym] = sscanf(line,"%d %s","C");
        z = mol_dbatom(atsym);
        idx = find(cr.ztyp == z);
        if (length(idx) > 1)
          idx = find(strcmp(cr.attyp,atsym));
        endif          
        if (idx > 0)
          cr.typ(i) = idx;
          line = substr(line,index(line,"=")+1);
          line = substr(line,index(line,"(")+1);
          cr.x(i,:) = sscanf(line,"%f %f %f");
          cr.x(i,:) = cr.x(i,:) * alat * rinv;
        else
          error(sprintf("Unknown atom number %d\n",i));
        endif
      endfor
    endif

    ## update with CELL_PARAMETERS and ATOMIC_POSITIONS
    ## (for crashed calculations)
    if(regexp(line,"^CELL_PARAMETERS"))
      [dum1 dum2 alat] = sscanf(line,"%s %s %s","C");
      alat = alat(1:length(alat)-1);
      if (!isempty(alat)) 
        alat = str2num(alat);
      else
        if (regexp(dum2,"bohr"))
          alat = 1;
        else
          error("fixme: cannot handle angstrom in CELL_PARAMETERS of the output")
        endif
      endif
      r = reshape(fscanf(fid,"%f"),3,3)';
      rinv = inv(r);
      nocrystalaxes = 1;
    endif
    if (regexp(line,"^ATOMIC_POSITIONS"))
      noatpos = 1;
      if (regexp(lower(line),"angstrom"))
        idocart = bohrtoans;
      endif
      cr.x = zeros(cr.nat,3);
      for i = 1:cr.nat
        line = fgetl(fid);
        [dum, cr.x(i,1), cr.x(i,2), cr.x(i,3)] = sscanf(line,"%s %f %f %f","C");
      endfor
    endif
    line = fgetl(fid);
  until (!ischar(line) && (line == -1))

  if (idocart > 0)
    rinv = inv(r);
    cr.x /= idocart;
    for i = 1:cr.nat
      cr.x(i,:) = cr.x(i,:) * rinv;
    endfor
  endif

  cr.r = r;
  cr.g = r * r';
  cr.a = sqrt(diag(cr.g))';
  cr.b(1) = acos(cr.g(2,3) / (cr.a(2)*cr.a(3)));
  cr.b(2) = acos(cr.g(1,3) / (cr.a(1)*cr.a(3)));
  cr.b(3) = acos(cr.g(1,2) / (cr.a(1)*cr.a(2)));
  cr.omega = det(cr.r);

  if (LOG > 0)
    printf("cr_read_espresso_out: Reading %s pwscf output\n", file);
    cr_popinfo(cr)
  endif

endfunction

