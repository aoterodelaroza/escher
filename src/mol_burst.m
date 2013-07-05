% Copyright (c) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
% Adapted from a tessel routine.
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

function smol = mol_burst(mol, bondfactor=1.20)
% function smol = mol_burst(mol, bondfactor=1.20)
%
% mol_burst - from a molecule (mol), extract a cell array of molecules
% containing each molecualr motif.
%
% Required input variables:
% mol: the molecule
%
% Optional input variables (all have default values):
% bondfactor: atoms are connected if they are at a distance less than
% bondfactor times the sum of covalent radii.
%
% Required output variables:
% {smol}: cell array of molecules containing the molecular motifs.
%

  ## distance matrix
  nat = length(mol.atnumber);
  lat = 1:nat;
  active = ones(1,nat);
  d = mol_distmatrix(mol);

  ## covalent radii
  rad = zeros(1,nat);
  for i = 1:nat
    [symb,atom] = mol_dbsymbol(mol.atnumber(i));
    rad(i) = atom.rcov;
  endfor

  ## build the fragments
  smol = cell();
  nout = 0;
  slat = lat(find(active));
  do
    idx = [slat(1)];
    active(slat(1)) = 0;

    do
      isnew = 0;
      for ii = 1:length(idx)
        i = idx(ii);
        for jj = 1:length(slat)
          j = slat(jj);
          if (!active(j))
            continue
          endif
          x = bondfactor*(rad(i)+rad(j));
          if (d(i,j) < x)
            active(j) = 0;
            idx = [idx j];
            isnew = 1;
          endif
        endfor
      endfor
    until (!isnew)

    slat = lat(find(active));
    nout += 1;
    smol{nout} = mol_getfragment(mol,idx);
  until (isempty(slat))

endfunction
