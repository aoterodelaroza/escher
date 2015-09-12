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

function rep = cr_unitcell(cr, addto="", i0=[0 0 0], i1=[0 0 0], radius=0.03, rgb=[255 0 0], tex="stick_default", axcolor="")
% function rep = cr_unitcell(cr, addto="", i0=[0 0 0], i1=[0 0 0], radius=0.05, rgb=[255 0 0], tex="stick_default", axcolor="")
%
% cr_unitcell - create sticks for the crystal unit cell.
%
% Required input variables:
% cr: structure with the crystal description.
%
% Optional input variables:
% addto: input representation on which the new representation is 
%        added.
% i0: initial lattice vectors. 
% i1: final lattice vectors. Represent all cells from i0 to i1.
%     Default: only the main cell (i0 = i1 = [0 0 0])
% radius: the stick radius. By default, 0.05.
% rgb: the stick color. By default, red. rgb can be given as three integer numbers 
%      (from 0 to 255) representing the rgb components. In addition, a fourth 
%      component (filter) and a fifth component (transmit) control the transparency
%      of the stick.
% tex: string identifier of the stick texture. This is interpreted (in subsequent calls
%      to the rep routines, i.e., not immediately) as the texture of the stick by calling 
%      the internal texture database.
% axcolor: if 1, plot only the x, y and z axes and use colors for them.
%

  bohr2ang = 0.52917720859;

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation_();
    if (isfield(cr,"name") && !isempty(cr.name))
      rep.name = cr.name;
    endif
  endif

  ## crystal to cartesian
  r = cr.r;
  g = cr.g;
  r *= bohr2ang;

  if (isfield(cr,"name"))
    name = strcat(cr.name,"_cell");
  else
    name = "cell";
  endif

  [rep itex] = rep_registertexture(rep,tex);

  if (isempty(axcolor)) 
    x0 = [
          0 0 0 1 0 0
          0 0 0 0 1 0
          0 0 0 0 0 1
          1 0 0 1 1 0
          1 0 0 1 0 1
          0 1 0 0 1 1
          0 1 0 1 1 0
          0 0 1 1 0 1
          0 0 1 0 1 1
          1 1 1 1 1 0
          1 1 1 1 0 1
          1 1 1 0 1 1
          ];
    for ix = i0(1):i1(1)
      for iy = i0(2):i1(2)
        for iz = i0(3):i1(3)
          for i = 1:size(x0,1)
            ## add the stick
            rep.nstick = rep.nstick + 1;
            rep.stick{rep.nstick} = stick();
            rep.stick{rep.nstick}.name = name;
            rep.stick{rep.nstick}.x0 = ([ix iy iz] + x0(i,1:3)) * r;
            rep.stick{rep.nstick}.x1 = ([ix iy iz] + x0(i,4:6)) * r;
            rep.stick{rep.nstick}.r = radius;
            rep.stick{rep.nstick}.rgb = rgb;
            rep.stick{rep.nstick}.tex = itex;
          end
        end
      end
    end
  else
    ## x
    rep.nstick = rep.nstick + 1;
    rep.stick{rep.nstick} = stick();
    rep.stick{rep.nstick}.name = name;
    rep.stick{rep.nstick}.x0 = [0 0 0] * r;
    rep.stick{rep.nstick}.x1 = [1 0 0] * r;
    rep.stick{rep.nstick}.r = radius;
    rep.stick{rep.nstick}.rgb = [255 0 0];
    rep.stick{rep.nstick}.tex = itex;
    ## y
    rep.nstick = rep.nstick + 1;
    rep.stick{rep.nstick} = stick();
    rep.stick{rep.nstick}.name = name;
    rep.stick{rep.nstick}.x0 = [0 0 0] * r;
    rep.stick{rep.nstick}.x1 = [0 1 0] * r;
    rep.stick{rep.nstick}.r = radius;
    rep.stick{rep.nstick}.rgb = [0 255 0];
    rep.stick{rep.nstick}.tex = itex;
    ## z
    rep.nstick = rep.nstick + 1;
    rep.stick{rep.nstick} = stick();
    rep.stick{rep.nstick}.name = name;
    rep.stick{rep.nstick}.x0 = [0 0 0] * r;
    rep.stick{rep.nstick}.x1 = [0 0 1] * r;
    rep.stick{rep.nstick}.r = radius;
    rep.stick{rep.nstick}.rgb = [0 0 255];
    rep.stick{rep.nstick}.tex = itex;
  endif

endfunction
