% Copyright (c) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function mol = mol_readcube (cubefile, LOG=0)
% function mol = mol_readcube (cubefile, LOG=0)
%
% mol_readcube - read in a gaussian cube file.
%
% Required input variables:
% {cubefile}: name of the cube file.
%
% Output:
% mol: molecule description. 
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2012

  bohrtoans = 0.52917720859;

  [fqub,msg] = fopen(cubefile,'r');
  if (fqub < 0 || ferror(fqub))
    disp(msg)
    error("mol_readcube: Could not open -- %s",cubefile);
  endif

  mol = molecule();

  ## title lines in the cube file
  mol.name = fgetl(fqub);
  title{2} = fgetl(fqub);

  ## read in number of atoms, origin, and grid
  line = fgetl(fqub);
  [tok,line] = strtok(line); mol.nat = str2num(tok);
  for i = 1:3
    line = fgetl(fqub);
  endfor

  ## read in atoms and atomic positions
  for i = 1:mol.nat
    line = fgetl(fqub);
    [tok,line] = strtok(line); mol.atnumber(i) = str2num(tok);
    mol.atname{i} = mol_dbsymbol(mol.atnumber(i));
    [tok,line] = strtok(line); 
    mol.atxyz(1:3,i) = str2num(line) * bohrtoans;
  endfor
  mol = mol_fillatmass(mol);

  fclose(fqub);

  if (LOG>0)
    printf('mol_readcube:\n');
    printf('File read: %s\n', cubefile);
    printf('Number of atoms: %d\n', mol.nat);
    printf('--Z-- ---coordinates---\n');
    for i = 1:mol.nat
      printf('%5d %5s ', mol.atnumber(i), mol.atname{i});
      printf('%15.9f', mol.atxyz(1:3,i));
      printf('\n');
    endfor
  endif

endfunction
