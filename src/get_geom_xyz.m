#! /usr/bin/octave -q
% Copyright (C) 2010 Victor Lua~na and Alberto Otero-de-la-Roza
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

function geom = get_geom_xyz (filexyz,LOG=1)
% function geom = get_geom_xyz(filexyz,LOG=1)
%
% get_geom_xyz - Produces the significative geometry for the molecule
% described in a xyz file.
%
% Required input variables:
% filexyz: name of the input data file in xyz format.
%
% % Optional input variables (all have default values):
% {LOG}: current level of printing. 0: no printing.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: October 2012

err = mol_dbstart(0);

mol = mol_readxyz(filexyz,:,0);

# Check the internal geometry?
mol = mol_internalgeometry(mol);

endfunction

## Use it as a script
if (!exist("argn"))
   if (nargin > 0)
      args = argv();
      for i = 1 : nargin
         get_geom_xyz(args{i});
      endfor
   else
      printf('Use as: get_xyz_geom_c6x6_yz3 file(s)\n');
   endif
endif
