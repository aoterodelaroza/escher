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
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
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
              1,op_dist(mol.atxyz(:,1),mol.atxyz(:,2)));
    elseif (i == 3)
      fprintf(fid,"%s %d %.10f %d %.10f\n",mol.atname{i},...
              2,op_dist(mol.atxyz(:,2),mol.atxyz(:,3)),...
              1,op_angle(mol.atxyz(:,1),mol.atxyz(:,2),mol.atxyz(:,3)));
    else 
      fprintf(fid,"%s %d %.10f %d %.10f %d %.10f\n",mol.atname{i},...
              i-1,op_dist(mol.atxyz(:,i-1),mol.atxyz(:,i)),...
              i-2,op_angle(mol.atxyz(:,i-2),mol.atxyz(:,i-1),mol.atxyz(:,i)),...
              i-3,op_dihedral(mol.atxyz(:,i-3),mol.atxyz(:,i-2),mol.atxyz(:,i-1),mol.atxyz(:,i)));
    endif
  endfor
  if (!isempty(file))
    fclose(fid);
  endif

endfunction

