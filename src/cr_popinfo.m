% Copyright (c) 2012 Alberto Otero-de-la-Roza and Victor Lua~na
% Adapted from an AOR (2011) routine.
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

function cr_popinfo(cr,file="")
% function cr_popinfo(cr,file="")
%
% cr_popinfo - output some information about the crystal cr.
%
% Required input variables:
% {cr}: crystal (structure)**2.
%
% Optional input variables (all have default values):
% {file}: output filename. If no file is present, stdout.
%

  printf("cell lengths (bohr): %.6f %.6f %.6f\n", cr.a(1:3));
  printf("cell angles (degs.): %.6f %.6f %.6f\n", cr.b(1:3) * 180/pi);
  printf("\n");
  printf("Atoms in the unit cell (%d):\n", cr.nat);
  printf("---i-- ----x----- ----y----- ----z----- -atom- --Z--\n");
  for i = 1:cr.nat
    printf("%6d %10.6f %10.6f %10.6f %6s %5d\n", i, cr.x(i,:), cr.attyp{cr.typ(i)}, cr.ztyp(cr.typ(i)));
  endfor

endfunction
