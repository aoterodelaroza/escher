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

function rep = rep_setbgcolor(repi,rgb,LOG=0);
% function rep = rep_setbgcolor(repi,rgb,LOG=0);
%
% rep_setbgcolor - set the background color of the representation
%
% Required input variables:
% repi: input representation.
% rgb: background color (three integer numbers).
%
% Optional input variables (all have default values):
% {LOG}: verbose level (0=silent,1=verbose).
%
% Required output variables:
% rep: output representation.

  rep = repi;

  rep.bgcolor = rgb(1:3);

  if (LOG>0) 
    printf("Background color set: <%d,%d,%d>",rgb(1:3));
  endif

endfunction
