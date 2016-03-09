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

function [cm,mass] = mol_cmass(mol,noweight=0)
% function [cm,mass] = mol_cmass(mol)
%
% mol_cmass - determine the center of mass of a molecule. If noweight=1,
% make all the atoms weigh the same.
%
% Required input variables:
% mol: the molecule.
% noweight: make all the atomic weights equal to 1. This 
% is useful when computing the geometric barycenter of the molecule. 
%
% Required output variables:
% cm: column vector with the center of mass coordinates.
% mass: total mass of the molecule or fragment.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: December 2011

  if (!noweight)
    mass = sum(mol.atmass);
    cm = sum(mol.atxyz * mol.atmass',2) / mass;
  else
    mass = 1;
    cm = sum(mol.atxyz,2) / mol.nat;
  endif

endfunction
