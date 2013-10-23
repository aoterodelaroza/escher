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

function rep = rep_addlight(repi,pos,r="",shadowless=0,intensity=1,LOG=0);
% function rep = rep_addlight(repi,pos,r="",shadowless=0,intensity=1,LOG=0);
%
% rep_addlight - add a light to the graphical representation
%
% Required input variables:
% repi: input representation.
% pos: position of the light.
%
% Optional input variables (all have default values):
% r: modelview transformation matrix.
% {LOG}: verbose level (0=silent,1=verbose).
%
% Required output variables:
% rep: output representation.

  rep = repi;
  if (!isfield(rep,"nlight"))
    rep.nlight = 0;
    rep.light = cell();
  endif

  rep.nlight = rep.nlight + 1;
  rep.light{rep.nlight}.x = pos;
  rep.light{rep.nlight}.color = [255 255 255];
  rep.light{rep.nlight}.shadowless = shadowless;
  rep.light{rep.nlight}.intensity = intensity;
  if (!isempty(r))
    rep.light{rep.nlight}.matrix = r;
  endif

  if (LOG>0) 
    printf("Added light number %d at <%.5f,%.5f,%.5f> angstrom, white\n",rep.nlight,pos)
  endif

endfunction
