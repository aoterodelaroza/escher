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

function  err = mol_writexyz (mol, filename="none", mode="axyz")
% function err = mol_writexyz (mol, filename="none", mode="axyz")
%
% mol_writexyz - write in the data of a molecule to a xyz file. Several
%               similar xyz formats are supported.
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
% {mode = "axyz"}: Several xyz formats are used by different programs.
%       Available possibilities are:
%       "axyz" ----> atomic symbol + (x,y,z) coordinates    
%       "nxyz" ----> atomic number + (x,y,z) coordinates    
%       "anxyz" ---> atomic symbol + atomic number + (x,y,z) coordinates    
%       "naxyz" ---> atomic number + atomic symbol + (x,y,z) coordinates    
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
%
% Required output variables:
% err: error code. 0 means no error.
%
% Format of the xyz data file:
% mode: Several xyz formats are used by different programs.
%       This routine is designed to guess the correct format.
%       Available possibilities are:
%       "axyz" ----> atomic symbol + (x,y,z) coordinates    
%       "nxyz" ----> atomic number + (x,y,z) coordinates    
%       "anxyz" ---> atomic symbol + atomic number + (x,y,z) coordinates    
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: June 2011

  angtobohr = 1.88972613288564;
  bohrtoans = 0.52917720859;

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

  nat = mol.nat;
  fprintf(fid," %d\n", nat);
  if (isfield(mol,"name") && !isempty(mol.name))
    fprintf(fid,"# %s\n",mol.name);
  else
    fprintf(fid,"# File title\n");
  endif
  for i = 1 : nat
    iz = mol_dbatom(mol.atname{i});
    sy = mol_dbsymbol(iz);
    if (strcmpi(mode, "axyz"))
      fprintf(fid," %-2s  ", sy);
      fprintf(fid," %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
      fprintf(fid,"\n");
    elseif (strcmpi(mode, "anxyz"))
      fprintf(fid," %-2s  ", sy);
      fprintf(fid," %-2d  ", iz);
      fprintf(fid," %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
      fprintf(fid,"\n");
    elseif (strcmpi(mode, "naxyz"))
      fprintf(fid," %-2d  ", iz);
      fprintf(fid," %-2s  ", sy);
      fprintf(fid," %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
      fprintf(fid,"\n");
    elseif (strcmpi(mode, "nxyz"))
      fprintf(fid," %-2d  ", iz);
      fprintf(fid," %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
      fprintf(fid,"\n");
    elseif (strcmpi(mode, "xyz"))
      fprintf(fid," %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
      fprintf(fid,"\n");
    else
      error("mol_writexyz", "Unknown write format!");
    endif
###sy = mol.atname{i};
###iz = mol_dbatom(sy);
###if (strcmpi(mode(1),"a"))
###   fprintf(fid," %-2s  ", sy);
###elseif (strcmpi(mode(1),"n"))
###   fprintf(fid," %-2d  ", iz);
###endif
###if (strcmpi(mode(2),"a"))
###   fprintf(fid," %2s  ", sy);
###elseif (strcmpi(mode(2),"n"))
###   fprintf(fid," %-2d  ", iz);
###endif
###fprintf(fid," %+15.9f %+15.9f %+15.9f", mol.atxyz(1:3,i));
###fprintf(fid,"\n");
  endfor

  if (!strcmpi(filename,"none"))
     fclose(fid);
  endif
  err = 0;

endfunction
