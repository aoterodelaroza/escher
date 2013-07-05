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

function rep = mol_stick(mol, addto="", s1=".+", s2=".+", dist=[-1 1.15], strict=1, radius=0.05, rgb=[255 0 0], tex="stick_default", LOG=0)
% function rep = mol_stick(mol, addto="", s1="", s2="", dist=0, strict=1, radius=0.05, rgb=[255 0 0], tex="stick_default", LOG=0)
%
% mol_stick - create sticks for a pair of atomic types given by their symbol.
% Optionally, use a distance criterion to build all the covalent bonds
% automatically.
%
% Required input variables:
% mol: structure with the molecular description.
%
% Optional input variables (all have default values):
% addto: input representation on which the new representation is 
%        added.
% s1: atomic symbol for the first atom type. By default, consider all
%     types. If strict = 1 (see below), a regular expression can be used. 
% s2: atomic symbol for the second atom type. By default, consider all
%     types. If strict = 1 (see below), a regular expression can be used.
% dist: minimum and maximum atomic distance to place a stick, in the form
%       dist = [dmin dmax]. If the form [-1 f] is used, then place a stick 
%       if the atomic distance is less than f times the sum of covalent radii.
% strict: if false (0), then s1 and s2 are converted to atomic numbers and 
%         sticks are represented for all atoms in mol with the same Z. Therefore,
%         s1 and s2 in this case should be plain strings.
%         If strict is true (1), then only atoms matching those names are used. 
%         Octave regular expressions (compared to the atomic names using regexp)
%         can be used in this case.
%         If strict is 2, use s1 and s2 as strings (compare with strcmp),
%         not regular expressions.
% radius: the stick radius. By default, 0.05.
% rgb: the stick color. By default, red. rgb can be given as three integer numbers 
%      (from 0 to 255) representing the rgb components. In addition, a fourth 
%      component (filter) and a fifth component (transmit) control the transparency
%      of the stick.
% tex: string identifier of the stick texture. This is interpreted (in subsequent calls
%      to the rep routines, i.e., not immediately) as the texture of the stick by calling 
%      the internal texture database.
% {LOG}: verbose level (0=silent,1=verbose).
%

  ## number of atoms
  nat = length(mol.atname);

  ## distance matrix
  nat = length(mol.atnumber);
  d = mol_distmatrix(mol);

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
    if (isfield(mol,"name") && !isempty(mol.name))
      rep.name = mol.name;
    endif
  endif
  if (!isfield(rep,"nstick"))
    rep.nstick = 0;
  endif
  if ((!isfield(rep,"stick")) || rep.nstick == 0)
    rep.stick = cell();
  endif

  ## covalent radii
  nat = length(mol.atnumber);
  if (dist(1) < 0)
    rcov = zeros(1,nat);
    for i = 1:nat
      [sym atom] = mol_dbsymbol(mol.atnumber(i));
      rcov(i) = atom.rcov;
    endfor
  endif

  ## connectivity matrix
  z1 = z2 = -1;
  if (!strict)
    z1 = mol_dbatom(s1,LOG);
    z2 = mol_dbatom(s2,LOG);
  endif
  for i = 1:nat
    for j = 1:i-1
      if (z1 > 0 && z2 > 0 || strict)
        ## are the atoms correct?
        if (strict == 2)
          doit = strcmp(mol.atname{i},s1) && strcmp(mol.atname{j},s2);
          doit = doit || (strcmp(mol.atname{i},s2) && strcmp(mol.atname{j},s1));
        elseif (strict == 1) 
          doit = regexp(mol.atname{i},s1) && regexp(mol.atname{j},s2);
          doit = doit || (regexp(mol.atname{i},s2) && regexp(mol.atname{j},s1));
        else
          zi = mol.atnumber(i);
          zj = mol.atnumber(j);
          doit = (zi == z1) && (zj == z2); 
          doit = doit || ((zi == z2) && (zj == z1)); 
        endif
        if (!doit) 
          continue
        endif
      endif

      ## and the distance?
      if (dist(1) < 0)
        doit = (d(i,j) <= dist(2) * (rcov(i) + rcov(j)));
      else
        doit = (d(i,j) >= dist(1)) && (d(i,j) <= dist(2));
      endif
      if (!doit) 
        continue
      endif

      ## add the stick
      rep.nstick = rep.nstick + 1;
      rep.stick{rep.nstick}.name = strcat(mol.atname{i},"_",mol.atname{j});
      rep.stick{rep.nstick}.x0 = mol.atxyz(:,i)';
      rep.stick{rep.nstick}.x1 = mol.atxyz(:,j)';
      rep.stick{rep.nstick}.r = radius;
      rep.stick{rep.nstick}.rgb = rgb;
      rep.stick{rep.nstick}.tex = tex;
    endfor
  endfor

endfunction
