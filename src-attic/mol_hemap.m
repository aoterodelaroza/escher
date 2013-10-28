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

function [map] = mol_hemap(file_xyz, file_cube="test.cube", LOG=1)
% function [map] = mol_hemap(file_xyz, file_cube="test.cube", LOG=1)
%
% mol_hemap - create the datafile for a He-map of the molecule in the
% xyz file.
%
% Required input variables:
% file_xyz: input molecule in the xyz format.
%
% Optional input variables (all have default values):
% {file_cube="test.cube"}: output map data to this file.
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%

mol_dbstart();

mol = mol_readxyz(file_xyz);

# Get the bounding box for the molecule:
BB = [min(mol.atxyz'); max(mol.atxyz')];

# The working box adds some extra space outside BB:
extra_space = 3.0;
WB = [BB(1,:)-extra_space; BB(2,:)+extra_space];

# We want a uniform grid on all the three directions. The step is fixed
# by the number of points in the shortest direction of the WB:
lenWB = WB(2,:)-WB(1,:);
cenWB = (WB(1,:) + WB(2,:)) / 2;
minWB = min(lenWB);
npt = 21;
###step = 1e-3 * round(1e3 * minWB / (npt-1));
step = minWB / (npt-1);

# We want an odd number of points on each direction (even intervals):
nptWB = round(lenWB / step);
even = rem(nptWB + 1,2);
nptWB = nptWB + even;
if (LOG > 1)
   printf("Grid size: %dx%dx%d (%d points)\n", nptWB, prod(nptWB));
endif

# The grid is defined from the centerpoint:
origWB = cenWB - ((nptWB-1)/2) * step;
for i = 1:nptWB(1)
   x = cenWB(1) + (i - (nptWB(1)+1)/2) * step;
   for j = 1:nptWB(2)
      y = cenWB(2) + (j - (nptWB(2)+1)/2) * step;
      for k = 1:nptWB(3)
         z = cenWB(3) + (k - (nptWB(3)+1)/3) * step;
         map{i,j,k}.xyz = [x,y,z];
         if (LOG > 1)
            printf("%4d%4d%4d", i, j, k);
            printf("%10.5f%10.5f%10.5f", x, y, z);
            printf("\n");
         endif
      endfor
   endfor
endfor

endfunction
