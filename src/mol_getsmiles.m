% Copyright (C) 2013 Alberto Otero-de-la-Roza and Victor Lua~na
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

function mol = mol_getsmiles(smiles)
% function mol = mol_getsmiles(smiles)
%
% mol_getsmiles - Use open babel to convert a smiles code into a molecule.
%
% smiles: a string containing the smiles code.

   [s out] = system(sprintf("obabel -:'%s' -oxyz -h --gen3d 2>&1", smiles));
   if (s)
     error(sprintf("Error using command: %s",order))
   endif

   mol = molecule_();
   mol.name = smiles;
   atext = textscan(out,"%s");
   atext = atext{1};
   mol.nat = str2num(atext{1});

   n = 2;
   for i = 1:mol.nat
     mol.atname{i} = atext{++n};
     for j = 1:3
       mol.atxyz(j,i) = str2num(atext{++n});
     endfor

     [atnumber,atprop] = mol_dbatom(mol.atname{i});
     mol.atnumber(i) = atnumber;
     mol.atmass(i) = atprop.mass;
   endfor

endfunction
