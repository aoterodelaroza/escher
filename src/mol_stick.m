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

function rep = mol_stick(mol, addto="", s1=".+", s2=".+", dist=[-1 1.15], strict=0, radius=0.05, rgb=[-1 -1 -1], tex="stick_default",round=0)
% function rep = mol_stick(mol, addto="", s1="", s2="", dist=0, strict=0, radius=0.05, rgb=[-1 -1 -1], tex="stick_default",round=0)
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
%         s1 and s2 in this case should be plain strings (except for the default .+).
%         If strict is true (1), then only atoms matching those names are used. 
%         Octave regular expressions (compared to the atomic names using regexp)
%         can be used in this case.
%         If strict is 2, use s1 and s2 as strings (compare with strcmp),
%         not regular expressions.
% radius: the stick radius. By default, 0.05.
% rgb: the stick color. By default, half and half. rgb can be given as three integer numbers 
%      (from 0 to 255) representing the rgb components. In addition, a fourth 
%      component (filter) and a fifth component (transmit) control the transparency
%      of the stick.
%      If any of the components is negative, use the internal color scheme to
%      color half the stick with one atom's color and half with the other.
% tex: string identifier of the stick texture. This is interpreted (in subsequent calls
%      to the rep routines, i.e., not immediately) as the texture of the stick by calling 
%      the internal texture database.
% round: use rounded cylinders in povray (requires shapes.inc).
%

  global atdb
  if (!exist("atdb","var") || isempty(atdb))
    err = mol_dbstart();
    if (err != 0)
      error("mol_stick: the atomic database does not start right!");
    endif
  endif

  ## distance matrix
  nat = mol.nat;
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

  ## covalent radii
  nat = mol.nat;
  if (dist(1) < 0)
    rcov = zeros(1,nat);
    for i = 1:nat
      [sym atom] = mol_dbsymbol(mol.atnumber(i));
      rcov(i) = atom.rcov;
    endfor
  endif

  isstick = zeros(nat);
  ## distance matrix to connectivity matrix
  if (dist(1) < 0)
    isstick = d <= dist(2) * (rcov' * ones(1,nat) + ones(nat,1) * rcov);
  else
    isstick = (d >= dist(1)) & (d <= dist(2));
  endif

  ## apply end-atom criteria
  z1 = z2 = -1;
  if (!strict)
    if (strcmp(s1,".+"))
      z1 = -1;
    else
      z1 = mol_dbatom(s1);
    endif
    if (strcmp(s1,".+"))
      z2 = -1;
    else
      z2 = mol_dbatom(s2);
    endif
  endif
  if (strict == 2) 
    r1 = r2 = zeros(1,nat);
    for i = 1:nat
      r1(i) = strcmp(mol.atname{i},s1);
      r2(i) = strcmp(mol.atname{i},s2);
    endfor
    isstick &= ((r1' * ones(1,nat)) & (ones(nat,1) * r2) | ...
                (r2' * ones(1,nat)) & (ones(nat,1) * r1));
  elseif (strict == 1)
    r1 = r2 = zeros(1,nat);
    for i = 1:nat
      r1(i) = !isempty(regexp(mol.atname{i},s1));
      r2(i) = !isempty(regexp(mol.atname{i},s2));
    endfor
    isstick &= ((r1' * ones(1,nat)) & (ones(nat,1) * r2) | ...
                (r2' * ones(1,nat)) & (ones(nat,1) * r1));
  elseif (z1 > 0 && z2 > 0)
    isstick &= ((mol.atnumber' * ones(1,nat) == z1) & (ones(nat,1) * mol.atnumber == z2) | ...
                (mol.atnumber' * ones(1,nat) == z2) & (ones(nat,1) * mol.atnumber == z1));
  endif

  dohalf = any(rgb < 0);

  ## add the sticks
  [inew,jnew] = find(isstick);
  nnew = length(inew);
  if (round && dohalf)
    ifac = 4;
  elseif (dohalf)
    ifac = 2;
  else
    ifac = 1;
  endif
  newstick = cell(1,ifac * nnew);
  [rep itex] = rep_registertexture(rep,tex);
  l = 0;
  for k = 1:nnew
    i = inew(k); j = jnew(k);
    l++;
    if (dohalf) 
      if (!round)
        xhalf = (mol.atxyz(:,i)' + mol.atxyz(:,j)') / 2;
        newstick{l} = stick();
        newstick{l}.name = strcat(mol.atname{i},"_",mol.atname{j},"_1");
        newstick{l}.x0 = mol.atxyz(:,i)';
        newstick{l}.x1 = xhalf;
        newstick{l}.r = radius;
        newstick{l}.rgb = fillrgb(atdb.color(1:3,mol.atnumber(i))');
        newstick{l}.tex = itex;
        newstick{l}.round = round;
        l++;
        newstick{l} = stick();
        newstick{l}.name = strcat(mol.atname{i},"_",mol.atname{j},"_2");
        newstick{l}.x0 = xhalf;
        newstick{l}.x1 = mol.atxyz(:,j)';
        newstick{l}.r = radius;
        newstick{l}.rgb = fillrgb(atdb.color(1:3,mol.atnumber(j))');
        newstick{l}.tex = itex;
        newstick{l}.round = round;
      else
        newstick{l} = stick();
        newstick{l}.name = strcat(mol.atname{i},"_",mol.atname{j},"_1");
        newstick{l}.x0 = mol.atxyz(:,i)';
        newstick{l}.x1 = (mol.atxyz(:,i)' + mol.atxyz(:,j)') / 2;
        newstick{l}.r = radius;
        newstick{l}.rgb = fillrgb(atdb.color(1:3,mol.atnumber(i))');
        newstick{l}.tex = itex;
        newstick{l}.round = round;
        l++;
        newstick{l} = stick();
        newstick{l}.name = strcat(mol.atname{i},"_",mol.atname{j},"_2");
        newstick{l}.x0 = (mol.atxyz(:,i)' + mol.atxyz(:,j)') / 2;
        newstick{l}.x1 = mol.atxyz(:,j)';
        newstick{l}.r = radius;
        newstick{l}.rgb = fillrgb(atdb.color(1:3,mol.atnumber(j))');
        newstick{l}.tex = itex;
        newstick{l}.round = round;
        l++;
        newstick{l} = stick();
        newstick{l}.name = strcat(mol.atname{i},"_",mol.atname{j},"_2");
        newstick{l}.x0 = 0.6 * mol.atxyz(:,i)' + 0.4 * mol.atxyz(:,j)';
        newstick{l}.x1 = (mol.atxyz(:,i)' + mol.atxyz(:,j)') / 2;
        newstick{l}.r = radius;
        newstick{l}.rgb = fillrgb(atdb.color(1:3,mol.atnumber(i))');
        newstick{l}.tex = itex;
        newstick{l}.round = 0;
        l++;
        newstick{l} = stick();
        newstick{l}.name = strcat(mol.atname{i},"_",mol.atname{j},"_2");
        newstick{l}.x0 = (mol.atxyz(:,i)' + mol.atxyz(:,j)') / 2;
        newstick{l}.x1 = 0.4 * mol.atxyz(:,i)' + 0.6 * mol.atxyz(:,j)';
        newstick{l}.r = radius;
        newstick{l}.rgb = fillrgb(atdb.color(1:3,mol.atnumber(j))');
        newstick{l}.tex = itex;
        newstick{l}.round = 0;
      endif
    else
      newstick{l} = stick();
      newstick{l}.name = strcat(mol.atname{i},"_",mol.atname{j});
      newstick{l}.x0 = mol.atxyz(:,i)';
      newstick{l}.x1 = mol.atxyz(:,j)';
      newstick{l}.r = radius;
      newstick{l}.rgb = rgb;
      newstick{l}.tex = itex;
      newstick{l}.round = round;
    endif
  endfor
  rep.nstick = rep.nstick + ifac * nnew;
  rep.stick = [rep.stick, newstick];

  ## register shapes.inc for loading if using rounded cylinders
  if (round)
    rep.load.shapes = 1;
  endif

endfunction
