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

function cr = crystal()
% function cr = crystal()
%
% crystal - create an empty crystal structure and initialize the number
% of atoms to zero.
%
% Output:
% {cr}: the empty crystal structure with all the fields defined.
%

  cr.name = "";
  cr.nat = 0;
  cr.ntyp = 0;
  cr.attyp = cell();
  cr.rvdwtyp = cr.c6typ = cr.zvaltyp = cr.qtyp = cr.ztyp = [];
  cr.typ = cr.x = [];
  cr.a = cr.b = zeros(1,3);
  cr.omega = 0;

endfunction
