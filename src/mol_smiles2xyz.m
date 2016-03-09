% Copyright (C) 2013 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [fxyz,mol] = mol_smiles2xyz (smiles="O=C(O)C(N)C Ala", code="ala", LOG=1)
% function [fxyz,mol] = mol_smiles2xyz (smiles, code, LOG=1)
%
% mol_smiles2xyz - Use open babel to convert a smiles code into a xyz
% file.
%
% Required input variables:
% {smiles}:
% {code}: root for the xyz filename:
%    <code>.xyz ------> cartesian coordinates.
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  debug information.
%
% Output arguments:
% {fxyz}: name of the file with the xyz coordinates.
%
% Optional output arguments:
% {mol}: molecular structure.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: Jan 2013

   bohrtoans = 0.52917720859;

   fxyz = sprintf("%s.xyz", code);
   order = sprintf("obabel -:'%s' -oxyz -h --gen3d -O %s", smiles, fxyz);
   system(order);

   if (nargout > 1)
      mol = mol_readxyz(fxyz);
   endif

   if (LOG > 0)
      printf ("smiles code: %s\n", smiles);
      printf ("coordinates in: %s\n", fxyz);
   endif

endfunction
