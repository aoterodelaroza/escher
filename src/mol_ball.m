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

function rep = mol_ball(mol, addto="", symb=".+", strict=1, radius=-0.6, rgb=[-1 -1 -1], tex="ball_default", LOG=0)
% function rep = mol_ball(mol, addto="", symb=".+", strict=1, radius=-0.6, rgb=[-1 -1 -1], tex="ball_default", LOG=0)
%
% mol_ball - create balls for an atomic type given by its symbol.
%
% Required input variables:
% mol: structure with the molecular description.
% symb: atomic symbol. By default, consider all atom types. If strict = 1, 
%       a regular expression can be used (see below).
%
% Optional input variables (all have default values):
% addto: input representation on which the new representation is 
%        added.
% strict: if false (0), then symb is converted to an atomic number and 
%         balls are represented for all atoms in mol with the same Z. If strict 
%         is true (1), then only atoms matching that name are represented. Octave
%         regular expressions (compared to the atomic names using regexp)
%         can be used in this case.
%         strict=2 is the same as 1, but use symb as a strings (compare with strcmp), 
%         not a regular expression.
% radius: the ball radius in angstrom. If a negative number is given,
%         the ball radius is taken as abs(radius) times the covalent radius.
% rgb: the ball color. By default (or if any of the rgb components is negative, use 
%      the internal color scheme. rgb can be given as three integer numbers (from 0 to 
%      255) representing the rgb components. In addition, a fourth component (filter)
%      and a fifth component (transmit) control the transparency of the ball.
% tex: string identifier of the ball texture. This is interpreted (in subsequent calls
%      to the rep routines, i.e., not immediately) as the texture of the ball by calling 
%      the internal texture database.
% {LOG}: verbose level (0=silent,1=verbose).
%

  ## number of atoms
  nat = length(mol.atname);

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
    if (isfield(mol,"name") && !isempty(mol.name))
      rep.name = mol.name;
    endif
  endif
  if (!isfield(rep,"nball"))
    rep.nball = 0;
  endif
  if (rep.nball == 0)
    rep.ball = cell();
  endif

  ## Create balls
  if (strict != 1 && strict != 2)
    zz = mol_dbatom(symb,LOG);
  endif
  for i = 1:nat
    if (strict == 2)
      doit = strcmp(mol.atname{i},symb);
    elseif (strict == 1)
      doit = regexp(mol.atname{i},symb);
    else
      doit = (zz == mol.atnumber(i));
    endif
    if (doit) 
      n = rep.nball = rep.nball+1;
      rep.ball{n} = struct();
      rep.ball{n}.x = mol.atxyz(:,i)';
      rep.ball{n}.name = mol.atname{i};
      [dum, atom] = mol_dbsymbol(mol.atnumber(i),LOG);
      if (radius>0) 
        rep.ball{n}.r = radius;
      else
        rep.ball{n}.r = atom.rcov * abs(radius); 
      endif
      if (all(rgb>=0))
        rep.ball{n}.rgb = rgb;
      else
        rep.ball{n}.rgb = [atom.color 0 0];
      endif
      rep.ball{n}.tex = tex;
    endif
  endfor

endfunction
