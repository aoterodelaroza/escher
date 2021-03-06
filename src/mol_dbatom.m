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

function [Z,atom] = mol_dbatom (symbol)
% function [Z,atom] = mol_dbatom (symbol)
%
% mol_dbatom - given an atomic symbol, return the atomic number and
% properties.
%
% Required input variables:
% symbol: name of the element to be identified.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: June 2011
% Modified: Jan 2012: (VLC) symbols can have a numerical label like C2 or C_2

global atdb

% Check if the database has been initialized:
if (!exist("atdb","var") || isempty(atdb))
   mol_dbstart();
endif

% Now get the atom and its properties:
[err,iz] = ismember(upper(symbol),atdb.symbols);
if (err == 0)
  ssymbol = substr(upper(symbol),1,2);
  [err,iz] = ismember(ssymbol,atdb.symbols);
endif
if (err == 0)
  ssymbol = substr(ssymbol,1,1);
  [err,iz] = ismember(ssymbol,atdb.symbols);
endif
if (err == 0)
  error(strcat("mol_dbatom: atomic symbol ->", symbol, "<- not found!"));
endif
Z = iz;
atom.number = iz;
atom.symbol = symbol;
atom.rcov   = atdb.rcov(iz);
atom.mass   = atdb.mass(iz);
atom.color  = atdb.color(iz);

endfunction
