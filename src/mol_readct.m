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

function [mol] = mol_readct(mol0,filename)
% function [mol] = mol_readct(mol0filename)
%
% mol_readct - read the connectivity of a molecule in chemdraw
% connection table format (.ct).  This routine reads in the
% connectivity information (mol.adjl), but NOT the positions.
% The order of the atoms has to be the same as in mol.
%
% Required input variables:
% mol0: input molecule.
% filename: name of the ct file.
%
% Output variables:
% mol: output molecule, with connectivity information.

  [fid,msg] = fopen(filename,"r");
  if (fid < 0 || ferror(fid)) 
    disp(msg)
    error("mol_readxyz: Could not open: %s",filename);
  endif

  mol = mol0;

  ## line 1: title
  line = fgetl(fid);

  ## line 2: number of atoms and number of bonds
  line = fgetl(fid);
  [nat, nb] = sscanf(line,"%s %s",'C');
  nat = str2num(nat);
  if (mol.nat != nat)
    error("Number of atoms in mol and ct file are inconsistent")
  endif
  nb = str2num(nb);

  ## skip the atomic positions
  for i = 1:nat
    line = fgetl(fid);
  endfor

  ## prepare the adjacency list - connect all nodes to themselves
  mol.adjl = zeros(mol.nat+nb,3);
  for i = 1:mol.nat
    mol.adjl(i,:) = [i i 1];
  endfor

  ## read the bonds into the connectivity matrix
  n = mol.nat;
  for i = 1:nb
    line = fgetl(fid);
    [i1, i2] = sscanf(line,"%s %s",'C');
    n++;
    mol.adjl(n,:) = [str2num(i1) str2num(i2) 1];
  endfor

  fclose(fid);

endfunction
