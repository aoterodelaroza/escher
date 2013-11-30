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

function mol_writezmat(mol,file="")
% function zmat = mol_writezmat(mol, file="")
%
% mol_writezmat - write a zmatrix using the cartesian coordinates in mol.
%                  The order is not changed.
%
% Required input variables:
% {mol}: the molecule, including the cartesian coordinates..
%
% Optional input variables (all have default values):
% {file}: write the z-matrix description to a file. Default: stdout.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: March 2012

  if (!isempty(file))
    fid = fopen(file,"w");
  else
    fid = stdout();
  endif

  for i = 1:mol.nat
    if (i == 1) 
      fprintf(fid,"%s\n",mol.atname{i});
    elseif (i == 2)
      fprintf(fid,"%s %d %.10f\n",mol.atname{i},...
              1,_dist(mol.atxyz(:,1),mol.atxyz(:,2)));
    elseif (i == 3)
      fprintf(fid,"%s %d %.10f %d %.10f\n",mol.atname{i},...
              2,_dist(mol.atxyz(:,2),mol.atxyz(:,3)),...
              1,_ang(mol.atxyz(:,1),mol.atxyz(:,2),mol.atxyz(:,3)));
    else 
      fprintf(fid,"%s %d %.10f %d %.10f %d %.10f\n",mol.atname{i},...
              i-1,_dist(mol.atxyz(:,i-1),mol.atxyz(:,i)),...
              i-2,_ang(mol.atxyz(:,i-2),mol.atxyz(:,i-1),mol.atxyz(:,i)),...
              i-3,_dih(mol.atxyz(:,i-3),mol.atxyz(:,i-2),mol.atxyz(:,i-1),mol.atxyz(:,i)));
    endif
  endfor
  if (!isempty(file))
    fclose(fid);
  endif

endfunction

## distance
function d = _dist(x0, x1)
   d = norm(x1 - x0);
endfunction

## angle between two points
function a = _ang(x0, x1, x2)
   n10 = x0-x1; n10 = n10/norm(n10);
   n12 = x2-x1; n12 = n12/norm(n12);
   a = acos(dot(n10,n12)) * 180/pi;
endfunction

## dihedral of 3 points
function d = _dih(x0, x1, x2, x3)
   b1 = x1-x0; b2 = x2-x1; b3 = x3-x2;
   b12 = cross(b1, b2); b23 = cross(b2, b3);
   d = -atan2(norm(b2) * dot(b1, b23), dot(b12, b23)) * (180/pi);
endfunction

