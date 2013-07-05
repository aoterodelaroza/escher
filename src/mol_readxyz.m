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

function [mol] = mol_readxyz (filename, mode='axyz', LOG=0)
% function [mol] = mol_readxyz (filename, mode='axyz', LOG=0)
%
% mol_readxyz - read in the data of a molecule from a xyz file. Several
%               similar xyz formats are supported.
%
% Required input variables:
% filename: name of the data file.
%
% Optional input variables (all have default values):
% {mode = "axyz"}: Several xyz formats are used by different programs.
%       This routine is designed to use a format by default but
%       try to guess the correct format. Available possibilities are:
%       "axyz" ----> atomic symbol + (x,y,z) coordinates    
%       "nxyz" ----> atomic number + (x,y,z) coordinates    
%       "anxyz" ---> atomic symbol + atomic number + (x,y,z) coordinates    
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  number of points read in, volume and energy range.
%            LOG = 2  like 1 plus a complete list of the points read in.
%
% Format of the xyz data file:
% mode: Several xyz formats are used by different programs.
%       This routine is designed to guess the correct format.
%       Available possibilities are:
%       "axyz" ----> atomic symbol + (x,y,z) coordinates    
%       "nxyz" ----> atomic number + (x,y,z) coordinates    
%       "anxyz" ---> atomic symbol + atomic number + (x,y,z) coordinates    
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

  if (LOG>0)
    tic;
    printf("mol_readxyz: Reading input file -> %s\n", filename);
  endif
  angtobohr = 1.88972613288564;
  bohrtoans = 0.52917720859;

  ## Datafile can contain the units for some critical variables. Currently:
  ## * columns -----> column order of data. Default: "axyz".
  ## * length ------> default unit: angstrom
  ## Those keywords will appear in comments in the head to the datafile.
  ## keywords
  Keyw={"columns","length"};

  [fid,msg] = fopen(filename,"r");
  if (fid < 0 || ferror(fid)) 
    disp(msg)
    error("mol_readxyz: Could not find -- %s",filename);
  endif

  natoms = 0;
  mol.name = 'unknown';
  mol.atname = {};
  mol.atnumber = [];
  mol.atxyz = [];

  ## Read in data file:
  columns = "axyz";
  atoms_read = -1;
  nl = 1;
  while (!feof(fid))
    line = fgetl(fid);

    ## headers and keywords
    [first,count] = sscanf(line," %1c",1);
    if (first == "\0" || first == " " || count == 0) 
      continue
    endif
    if (first == "#")
      # interpret the header as a sequence of lines:
      #   #  keyword  val
      # with keyword necessarily one of the keyw cell array.
      if (nl == 1) 
        # remove # chars and trailing blanks and tabs
        while (line(1) == "#" || line(1) == " " || line(1) == "\t")
          line = line(2:end);
        endwhile    
        [first,sec] = sscanf(line,"%s %s","C");
        if (isa(first,"char") && isa(sec,"char"))
          first = tolower(first); sec = tolower(sec);
          if(ismember(first,Keyw))
            eval(strcat(first,"='",sec,"' ;"));
          endif
        endif
      endif
      continue    # Skip to next line
    endif

    ## read the line and advance
    #if (columns == "axyz")
    if (true)
       if (LOG > 1)
          printf("DBG(mol_readxyz): line -> %s\n", line);
       endif
       [g{1:5},count] = sscanf(line,"%s %s %s %s %s",'C');
       if (count==1 && nl==1)
          atoms_read = cellfun('str2num',g(1));
          mol.name = fgetl(fid);
          nl = 3;
          continue
       elseif (count==4)
          xyz = cellfun('str2num',g(2:4));
          [sv,st] = cellfun('str2double',g(1));
          if (st == 0)
             columns = "nxyz";
             atnumber = sv;
             [atsymbol,atprop] = mol_dbsymbol(atnumber);
          else
             columns = "axyz";
             atsymbol = cell2mat(g(1));
             [atnumber,atprop] = mol_dbatom(atsymbol);
          endif
       elseif (count==5)
          xyz = cellfun('str2num',g(3:5));
          [sv,st] = cellfun('str2double',g(1));
          if (st == 0)
             columns = "naxyz";
             atnumber = sv;
             atsymbol = cell2mat(g(2));
          else
             columns = "anxyz";
             [atnumber,st] = cellfun('str2double',g(2));
             atsymbol = cell2mat(g(1));
             if (st != 0)
                error(['mol_readxyz: error reading line ', line]);
             endif
          endif
          [symb,atprop] = mol_dbsymbol(atnumber);
       else
          error(['mol_readxyz: error reading line ', line]);
       endif
       natoms = natoms + 1;
       mol.atname{natoms} = atsymbol;
       mol.atnumber(natoms) = atnumber;
       mol.atmass(natoms) = atprop.mass;
       %%%mol.atxyz{natoms} = xyz';
       mol.atxyz(1:3,natoms) = xyz';
    else
       error("mol_readxyz: unknown order of data columns!");
    endif
    nl++;
  endwhile
  fclose(fid);

  if (LOG>0)
    printf("Data file format: %s\n", columns);
    printf("Number of atoms: %d\n", natoms);
    printf("number, symbol, at_number, xyz coordinates:\n");
    for i = 1:natoms
       printf("%4d %2s %3d", i, mol.atname{i}, mol.atnumber(i));
       %%%printf(" %12.6f %12.6f %12.6f", cell2mat(mol.atxyz(i)));
       printf(" %12.6f %12.6f %12.6f", mol.atxyz(1:3,i));
       printf("\n");
    endfor
    toc;
  endif

endfunction
