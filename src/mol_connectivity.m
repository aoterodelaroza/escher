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

function mol = mol_connectivity(mol0)
% function mol = mol_connectivity(mol0)
%
% mol_connectivity - generate the edge list representation of a molecule mol
% and save it to the mol.adjl field. This routine requires openbabel (obabel)
% to work.
%
% Required input variables:
% mol: input molecule.
%
% Output variables:
% mol: output molecule. mol.adjl contains the adjacency list.

  ## make sure openbabel is in the path
  mol = mol0;
  [s out] = system("which obabel");
  if (s)
    error("mol_connectivity requires openbabel in your PATH to work (obabel).")
  endif

  ## prepare temporary files
  root = tempname;
  fxyz = sprintf("%s.xyz",root);
  fct = sprintf("%s.ct",root);
  mol_writexyz(mol,fxyz);

  ## run openbabel
  [s out] = system(sprintf("obabel -ixyz %s -oct -O %s 2>&1 > /dev/null",fxyz,fct));
  if (s)
    error(sprintf("Error executing command: obabel -ixyz %s -oct -O %s",fxyz,fct))
  endif

  ## read in the connectivity information
  mol = mol_readct(mol,fct);

  ## delete the temporary files
  [s out] = unlink(fxyz);
  if (s)
    error(sprintf("Error deleting file: %s",fxyz))
  endif
  [s out] = unlink(fct);
  if (s)
    error(sprintf("Error deleting file: %s",fct))
  endif

endfunction
