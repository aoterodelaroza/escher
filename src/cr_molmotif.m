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

function [mol mask] = cr_molmotif(cr, ilist=[], bondfactor=1.20, LOG=0)
% function [mol mask] = cr_molmotif(cr, ilist=[], bondfactor=1.20, LOG=1)
%
% cr_molmotif - return one or several molecular motifs where, at least,
% one atom is within the main unit cell.
%
% Required input variables:
% {cr}: struct containing the crystal.
%
% Optional input variables (all have default values):
% {ilist[]}: list of atoms that will be used to determine molecular motifs.
%            Default: all the atoms in the unit cell.
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {mol}: description of the molecular motifs, including the cartesian
%      coordinates of the atoms.
% {mask}: cell array containing the mask that generates the mol. The
%      mask contains the atom indices and lattice translations.
%

global atdb

%% initialize
bohr2angstrom = 0.52917720859;
if (!exist("atdb","var") || isempty(atdb))
   err = mol_dbstart(LOG);
   if (err != 0)
      error("mol_dbatom: the atomic database does not start right!");
   endif
endif

%% Add the ilist vectors
mol = molecule();
mask = crmask();
if (isempty(ilist))
   ilist = 1:cr.nat;
endif

%% build the initial molecule
for i = ilist
    mol = mol_addatom(mol,cr.attyp{cr.typ(i)},cr.x(i,:) * cr.r * bohr2angstrom);
    mask.nat += 1;
    mask.l(mask.nat,:) = [0 0 0];
    mask.i(mask.nat) = i;
endfor

%% build the big molecule
mol0 = molecule();
mask0 = crmask();
for j = [0 -1 1]
  for k = [0 -1 1]
    for l = [0 -1 1]
      for i = 1:cr.nat
        mol0 = mol_addatom(mol0,cr.attyp{cr.typ(i)},(cr.x(i,:)+[j k l]) * cr.r * bohr2angstrom);
        mask0.nat += 1;
        mask0.l(mask0.nat,:) = [j k l];
        mask0.i(mask0.nat) = i;
      endfor
    endfor
  endfor
endfor

%% mark the initial list as done
mlist = zeros(1,mol0.nat);
mlist(ilist) = 1;

%% calculate the distance matrix
d = mol_distmatrix(mol0);

%% iteratively add to the molecule

do
  new = 0;
  newmlist = mlist;
  for i = find(!mlist)
    for j = find(mlist)
      if (d(i,j) <= bondfactor*atdb.rcov(mol0.atnumber(i))+atdb.rcov(mol0.atnumber(j)))
        newmlist(i) = 1;
      endif
    endfor
  endfor
  for i = find(newmlist - mlist)
    mol = mol_addatom(mol,mol0.atname{i},mol0.atxyz(:,i));
    mask.nat += 1;
    mask.l(mask.nat,:) = mask0.l(i,:);
    mask.i(mask.nat) = mask0.i(i);
  endfor
  if (any(newmlist - mlist))
    new = 1;
  endif
  mlist = newmlist;
until (!new)

endfunction
