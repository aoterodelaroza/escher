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

function limits = mol_limits (mol, mode='normal', LOG=0)
% function limits = mol_limits (mol, mode='normal', LOG=0)
%
% mol_limits - Determine the geometric limits of the molecule.
%
% Required input variables:
% mol: description of the molecule.
%
% Optional input variables (all have default values):
% {mode}: method for determining the limits. Currently coded:
%   * normal -> boundary box parallel to the xyz axes.
%
% Required output variables:
% limits: boundary box. The information depends on the method:
%   * normal -> limits.min(1:3) and limits.max(1:3) give the extension
%               along the xyz axes.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: January 2012

if (isequal(method,'normal'))
   limits.min = min(mol.atxyz')';
   limits.max = max(mol.atxyz')';
   if (LOG > 0)
      printf("mol_limits: bounding box of %s\n", mol.name);
      printf("min(x,y,z): %.6f %.6f %.6f\n", limits.min);
      printf("max(x,y,z): %.6f %.6f %.6f\n", limits.max);
      printf("\n");
   endif
else
   error('mol_limits: required method unknown!');
endir

endfunction
