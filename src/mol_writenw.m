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

function err = mol_writenw (mol, filename="none", LOG=0)
% function err = mol_writenw (mol, filename="none", LOG=0)
%
% mol_writenw - write in the data of a molecule to a nw file (input for
%               NWChem calculation).
%
% Required input variables:
% mol: structure with the molecular description. The format is:
%       mol.name --> name of the molecule.
%       mol.atname --> {1:M} cell array with the symbols of the atoms
%                       (M is the number of atoms in the molecule).
%       mol.atxyz --> Mx3 matrix with the cartesian coordinates of the atoms.
%
% Optional input variables (all have default values):
% {filename="none"): name of the output data file. By defult, the standard
%                    output is used.
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
%
% Required output variables:
% err: error code. 0 means no error.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: November 2011

  angtobohr = 1.88972613288564;
  bohrtoans = 0.52917720859;

  if (LOG > 0)
     mol
  endif

  err = 1;
  if (strcmpi(filename,"none"))
     fid = stdout();
  else
     [fid,msg] = fopen(filename,"w+");
     if (fid < 0 || ferror(fid)) 
       disp(msg)
       error("mol_readxyz: Could not find -- %s",filename);
     endif
  endif

  atcount = zeros(1,110);
  for i = 1 : mol.nat
     [ZZ,at] = mol_dbatom(mol.atname{i});
     atcount(ZZ)++;
     atsymbol{ZZ} = mol.atname{i};
  endfor
  # After "Lw" there can be a number of ghost and dummy atoms:
  [ZZ,at] = mol_dbatom("Lw");
  maxcount = ZZ;

  fprintf(fid, "start wkfiles\n");
  ###fprintf(fid, "memory total 8192 mb\n");
  ###fprintf(fid, "memory total 32000 mb\n");
  fprintf(fid, "memory total 64000 mb\n");
  #fprintf(fid, 'title "(VLC00) %s (mp2/6-311G**) PES exploration" \n', mol.name);
  fprintf(fid, 'title "(VLC00) %s (mp2/6-311G**) PES exploration 2" \n', mol.name);
  fprintf(fid, "charge 0\n");
  fprintf(fid, "geometry units angstroms print xyz\n");
  nat = mol.nat;
  for i = 1 : nat
     fprintf(fid,"   %-2s  ", mol.atname{i});
     fprintf(fid," %12.6f %12.6f %12.6f", mol.atxyz(1:3,i));
     fprintf(fid,"\n");
  endfor
  fprintf(fid, "end\n\n");
  fprintf(fid, "basis\n");
  # Be careful with ghost and dummy atoms
  for i = 1 : maxcount
     if (atcount(i) > 0)
        fprintf(fid,"   %-2s library 6-311G**\n", atsymbol{i});
     endif
  endfor
  fprintf(fid, "end\n\n");
  fprintf(fid, "mp2; tight; end\n");
  ###fprintf(fid, "driver; tight; maxiter 200; end\n");
  fprintf(fid,"\n");
  ###fprintf(fid, 'title "(VLC01: optim all) %s (mp2/6-311G**) PES exploration" \n', mol.name);
  fprintf(fid, 'title "(VLC01: energy) %s (mp2/6-311G**) PES exploration" \n', mol.name);
  ###fprintf(fid, "task mp2 optimize\n");
  fprintf(fid, "task mp2 energy\n");

  if (!strcmpi(filename,"none"))
     fclose(fid);
  endif
  err = 0;

endfunction
