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

function [symbol,atom] = mol_dbsymbol (Z,LOG=0)
% function [symbol,atom] = mol_dbsymbol (Z,LOG=0)
%
% mol_dbsymbol - given an atomic number, return the atomic name and
% properties.
%
% Required input variables:
% Z: atomic number of the element to be identified.
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  debug information.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: June 2011

   global atdb

%  Check if the database has been initialized:
   %%%if (!exist("atdb","var")  ||  atdb.defined!=true)
   %%%if (exist("atdb","var")!=1 | dbdefined!=1)
   if (!exist("atdb","var") || isempty(atdb))
     err = mol_dbstart(LOG);
     if (err != 0)
       error("mol_dbatom: the atomic database does not start right!");
     endif
   endif

   atom = atom();
%  Now get the atom and its properties:
   if (Z == 0)
     symbol = "Bq";
     atom.number = Z;
     atom.symbol = symbol;
     atom.mass   = 0;
     atom.rcov   = 0;
     atom.color  = [0 255 0];
     return
   elseif (Z < 1 || Z > length(atdb.rcov))
      error(["mol_dbsymbol: atomic number ->", int2str(Z), "<- not found!"]);
   endif
   symbol = cell2mat(atdb.symbols(Z));
   atom.number = Z;
   atom.symbol = symbol;
   atom.mass   = atdb.mass(Z);
   atom.rcov   = atdb.rcov(Z);
   atom.color  = atdb.color(1:3,Z)';

endfunction
