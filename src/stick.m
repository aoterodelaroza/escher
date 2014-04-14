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

function st = stick()
% function st = stick()
%
% stick - create an empty stick.
%
% Output:
% {st}: the empty stick.
%

  st.name = "";
  st.x0 = [0 0 0];
  st.x1 = [0 0 0];
  st.r = 0;
  st.rgb = [0 0 0 0 0];
  st.tex = 0;
  st.round = 0;

endfunction
