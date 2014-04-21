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

function run_povray(name,crop=1)
% function run_povray(name)
%
% run_povray - run povray on name.pov with the default parameters.
% If crop=1, then use imagemagick's "convert" to crop the margins.
%

  system(sprintf("povray -D -UV +I%s.pov +O%s.png +W2000 +H2000 +A",name,name));
  if (crop)
    system(sprintf("convert -trim -bordercolor White -border 0x0 +repage %s.png %s_crop.png",name,name));
    system(sprintf("mv %s_crop.png %s.png",name,name));
  endif

endfunction
