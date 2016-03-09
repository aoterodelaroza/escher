% Copyright (C) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [x0, nx, dx] = mol_cuberange (mol, qual=6, marg=1, cube=0, LOG=1)
% function [x0, nx, dx] = mol_cuberange (mol, quality=4, marg=1, cube=0, LOG=1)
%
% mol_cuberange - determines the range to use in cubegen.
% Notice: distances are assumed to be in angstrom.
%
% Required input variables:
%
% {mol}: mol to be calculated. Careful with the compatibility of the
% wavefunction. g09 has the bad habit to displace and reorient the input
% molecule except when the NoSymm option is being used.
%
% Optional input variables (all have default values):
% {qual}: density of points for the calculation grid.
%    Old version: levels 1 to 6 [3, 6, 12, 18, 24, 36] points/distance_units.
%    New version: quality is 3*qual points/distance_units.
%    (distance_units are usually bohr)
% {marg}: empty space or margin surrounding the finite molecule.
%    Default: 1 cell empty on each direction.
% {cube}: toggle that decides the use of a cube or just parallepipedic
%    grid. Default 0: non cube.
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  debug information.
%
% Required output variables:
% {x0}: origin for the grid.
% {nx}: number of points in the grid.
% {dx}: grid resolution.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: Dec 2012

   qual = fix(qual);
   if (qual >= 1)
      nx = 3*qual;
   else
      error("cuberange: unkown quality option - %d\n", quality);
   endif
   dx = 1.0/(nx-1);
   margin = [marg,marg,marg]';

   x0 = min(mol.atxyz(1:3,:)')' - margin;
   x1 = max(mol.atxyz(1:3,:)')' + margin;
   xr = x1 - x0;
   nx = fix(xr/dx+1);
   dx = xr./nx;

   if (cube > 0)
      x0min = min(x0);
      x1max = max(x1);
      xr = x1max - x0min;
      nxmax = max(nx);
      dxmin = xr / nxmax;
      x0 = ones(3,1) * x0min;
      x1 = ones(3,1) * x1max;
      nx = ones(3,1) * nxmax;
      dx = ones(3,1) * dxmin;
   endif

   if (LOG > 0)
      printf("%5d %12.6f %12.6f %12.6f\n", -6, x0);
      printf("%5d %12.6f %12.6f %12.6f\n", nx(1), dx(1), 0.0, 0.0);
      printf("%5d %12.6f %12.6f %12.6f\n", nx(2), 0.0, dx(2), 0.0);
      printf("%5d %12.6f %12.6f %12.6f\n", nx(3), 0.0, 0.0, dx(3));
   endif


endfunction
