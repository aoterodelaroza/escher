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

function [X0,X1] = mol_fragdim(mol, frag, LOG=1)
% function [X0,X1] = mol_fragdim(mol, frag, LOG=0)
%
% mol_fragdim- Returns the range (min and max coordinates)
%    of the entered fragment of the molecule.
%
% Required input variables:
% {mol}: The molecule to fragment.
% {frag}: fragment: Subset of the atoms in mol.
%
% Output:
% {X0,X1}: min and max values for the x, y, z, coordinates of the
%   atoms in the fragment of mol.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Aug 2015
% early attempt:
%  for i=1:3
%     X0(i) = min(mol.atxyz(i,frag));
%     X1(i) = max(mol.atxyz(i,frag));
%  endfor
  X0 = min(mol.atxyz(:,frag),[],2);
  X1 = max(mol.atxyz(:,frag),[],2);
  if (LOG>0)
     printf('mol_fragdim: Determine the orthohedron of the mol fragment\n');
     printf('num atoms in frag: %d\n', length(frag));
     printf('X0: %.5g %.5g %.5g\n', X0);
     printf('X1: %.5g %.5g %.5g\n', X1);
     if (LOG>1)
        for i = 1:length(frag)
           i1 = frag(i);
           printf('%3s ', mol.atname{1,i1});
           printf('%13.5f %13.5f %13.5f\n', mol.atxyz(1:3,i1));
        endfor
     endif
  endif

endfunction
