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

function [izmat izvar mol] = mol_sanezmatrix(mol0,dihlist)
% function [izmat izvar mol] = mol_sanezmatrix(mol0,dihlist)
%
% mol_sanezmatrix - write z-matrix in a format that allows rotation
% around the given list of dihedrals without changing the distances
% and angles in any of the atomic groups connected by those dihedrals.
% This routine is useful for systematic searches over conformers of
% small molecules.
%
% Required input variables:
% mol: input molecule.
% dihlist: [
%        i1 i2 i3 i4
%        i5 i6 i7 i8
%        ...]
%      Turn the dihedrals in this list.  These bonds have to be edges
%      in the connectivity graph. If they are not part of the graph,
%      you can precompute the connectivity graph before passing the
%      molecule and add them.
%
% Output variables:
% izmat: indices for the z-matrix, in a format understandable by mol_writezmat.
% izvar: variables for the z-matrix, in a format understandable by mol_writezmat.
%        The dihedrals around which the molecule rotates are included here.
%        The integer in the non-zero entries corresponds gives the line
%        in dihlist containing the corresponding dihedral.
% mol: output molecule, containing the atomic connectivity if it was
% not provided in input.
%

  ## calculate the adjacency list
  mol = mol0;
  if (!isfield(mol,"adjl") || isempty(mol.adjl))
    mol = mol_connectivity(mol);
  endif

  ## verify that all the bonds in the list are connected
  for i = 1:rows(dihlist)
    for j = 1:3
      idx = dihlist(i,j:j+1);
      if (!any(all(mol.adjl(:,1:2) == idx,2) | all(mol.adjl(:,1:2) == fliplr(idx),2)))
        error(sprintf("Could not find bond %d-%d in connectivity list",idx))
      endif
    endfor
  endfor

  ## check that the networks toolbox is present
  list = {"symmetrizeEdgeL","edgeL2adjL","adjL2adj","numConnComp","findConnComp","findConnCompI"};
  ex = cellfun("exist",list);
  if (!all(ex == 2 | ex == 5))
    disp("Missing graph theory routines!")
    disp("Please, install the octave-networks-toolbox package")
    disp("From: https://github.com/aeolianine/octave-networks-toolbox.git")
    error("Missing functions for mol_conformers.")
  endif

  ## save the dihlist
  dihlist0 = dihlist;

  ## make the graph an undirected graph
  ed = mol.adjl;
  eds = symmetrizeEdgeL(ed);

  ## adjacency matrix
  adjL = edgeL2adjL(eds);
  adj = adjL2adj(adjL);

  ## number of connected components
  nc = numConnComp(adj);
  if (nc > 1)
    error("FIXME: handle groups of disconnected molecules")
  endif

  ## find the pruned graph
  ped = [];
  for j = 1:rows(ed)
    found = 0;
    for i = 1:rows(dihlist)
      if (dihlist(i,2) == ed(j,1) && dihlist(i,3) == ed(j,2) || ...
                                     dihlist(i,3) == ed(j,1) && dihlist(i,2) == ed(j,2))
        found = 1;
      endif
    endfor
    if (!found)
      ped = [ped; ed(j,:)];
    endif
  endfor

  ## components of the pruned graph
  peds = symmetrizeEdgeL(ped);
  adjpL = edgeL2adjL(peds);
  adjp = adjL2adj(adjpL);
  ncp = numConnComp(adjp);

  if (ncp-1 != rows(dihlist))
    error("Can not rotate around one of the dihedrals -> cycle.")
  endif

  ## find the connected components; sort by decreasing size
  pcon = findConnComp(adjp);
  x = cellfun('length',pcon,1);
  [x idx] = sort(x,'descend');
  pcon = pcon(idx);

  ## check the sanity of the connected components
  if (length(pcon{1}) + length(pcon{2}) < 4)
    error("(fixme?) not enough atoms for a dihedral")
  endif

  ## Prepare the z-matrix and add the first connected component.
  dihsave = zeros(mol.nat,1);
  xcar = zeros(mol.nat,3);
  izmat = zeros(mol.nat,4);
  nat = 0;
  for i = 1:length(pcon{1})
    nat++;
    xcar(nat,:) = mol.atxyz(:,pcon{1}(i))';
    izmat(nat,:) = max([pcon{1}(i) nat-1 nat-2 nat-3],0);
  endfor

  while (!isempty(dihlist))
    ## find a dihlist entry contained in the already-built z-matrix
    rang = rows(dihlist);
    for i = 1:rows(dihlist)
      if (all(ismember(dihlist(i,1:2),izmat(1:nat,1)')))
        idx = dihlist(i,:);
        dihlist = dihlist([1:i-1, i+1:rang],:);
        break
      elseif (all(ismember(dihlist(i,3:4),izmat(1:nat,1)')))
        idx = fliplr(dihlist(i,:));
        dihlist = dihlist([1:i-1, i+1:rang],:);
        break
      endif
    endfor

    ## find the connected component for the new atoms
    con = findConnCompI(adjp,idx(3));

    ## add first connector atom 
    nat++;
    xcar(nat,:) = mol.atxyz(:,idx(3))';
    izmat(nat,:) = [idx(3) find(izmat(1:nat-1,1)==idx(2)) find(izmat(1:nat-1,1)==idx(1)) 0];
    con = setdiff(con,idx(3));

    ## for the dihedral: maybe it was saved from some other dihedral
    ii = find(izmat(:,1) == idx(1));
    if (dihsave(idx(3)) > 0)
      izmat(nat,4) = dihsave(idx(3));
    elseif (izmat(ii,2) > 0 && izmat(ii,2) != find(izmat(:,1) == idx(2)))
      ## maybe add the atom directly connected to idx(1)
      izmat(nat,4) = izmat(ii,2);
    elseif (nat > 3)
      ## find the connected component of idx(1)
      ncon = findConnCompI(adjp,idx(1));
      ifound = 0;
      for i = 1:length(ncon)
        if (ncon(i) != idx(1) && ncon(i) != idx(2))
          ifound = ncon(i);
          break
        endif
      endfor
      if (!ifound)
        error("Could not find connected dihedral for 1st connector; try reordering the molecule")
      endif
      izmat(nat,4) = find(izmat(1:nat-1,1)==ifound);
    endif

    ## add second connector atom, if it's part of this component
    if (ismember(idx(4),con))
      nat++;
      xcar(nat,:) = mol.atxyz(:,idx(4))';
      izmat(nat,:) = [idx(4) find(izmat(1:nat-1,1)==idx(3)) find(izmat(1:nat-1,1)==idx(2)) find(izmat(1:nat-1,1)==idx(1))];
      con = setdiff(con,idx(4));
    else
      dihsave(idx(4)) = find(izmat(1:nat-1,1)==idx(1));
    endif

    ## add the rest of the atoms
    for i = 1:length(con)
      nat++;
      xcar(nat,:) = mol.atxyz(:,con(i))';
      izmat(nat,:) = max([con(i) nat-1 nat-2 nat-3],0);
      if (i == 1)
        izmat(nat,4) = izmat(nat-2,2);
      endif
    endfor
  endwhile

  ## find the dihedrals 
  dihlist = dihlist0;
  izvar = zeros(size(izmat));
  for i = 1:rows(izmat)
    if (izmat(i,4) == 0)
      continue
    endif
    idx = [izmat(i,1) izmat(izmat(i,2),1) izmat(izmat(i,3),1) izmat(izmat(i,4),1)];
    found = 0;
    for j = 1:rows(dihlist)
      if (all(idx == dihlist(j,:)) || all(fliplr(idx) == dihlist(j,:)))
        found = j;
      endif
    endfor
    if (found > 0)
      izvar(i,4) = found;
    endif
  endfor

endfunction
