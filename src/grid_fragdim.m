% Copyright (c) 2015 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [X0,X1,gridX0,gridX1] = grid_fragdim(mol, grid, frag, LOG=0)
% function [X0,X1,gridX0,gridX1] = grid_fragdim(mol, grid, frag, LOG=0)
%
% grid_fragdim- Returns the range (min and max coordinates) of the entered fragment
%   of the mol, and the whole X0 and X1 range of the grid.
%
% Required input variables:
% {mol}: the molecule to fragment.
% {grid}: grid of the molecule to fragment.
% {frag}: fragment: Subset of the atoms in mol.
%
% Output:
% {X0,X1}: min and max values for the x, y, z, coordinates of the
%   atoms in the fragment of mol.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0. Print the fragment atom coordinates ig LOG>1.
%
% Example:
%   mol = mol_readcube('anthrazene-rho.cube);
%   rho = grid_readcube('anthrazene-rho.cube);
%   r1 = [1    2    3    4   13   14];  # the C atoms in the first ring
%   [X0,X1,gX0,gX1] = mol_fragdim(mol,rho,r1);
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: Aug 2015

  X0 = min(mol.atxyz(:,frag),[],2);
  X1 = max(mol.atxyz(:,frag),[],2);
  gridX0 = grid.x0';
  gridX1 =  (grid.x0 + grid.n.*diag(grid.dx)')';
  if (LOG>0)
     printf('mol_fragdim: Determine the orthohedron of the mol fragment\n');
     printf('num atoms in frag: %d\n', length(frag));
     printf('X0....: %.5g %.5g %.5g\n', X0);
     printf('X1....: %.5g %.5g %.5g\n', X1);
     printf('gridX0: %.5g %.5g %.5g\n', gridX0);
     printf('gridX1: %.5g %.5g %.5g\n', gridX1);
     if (LOG>1)
        for i = 1:length(frag)
           i1 = frag(i);
           printf('%3s %3d ', mol.atname{1,i1},i1);
           printf('%13.5f %13.5f %13.5f\n', mol.atxyz(1:3,i1));
        endfor
     endif
  endif

endfunction
