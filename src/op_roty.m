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

function m = op_roty (angle)
% function m = op_roty (angle)
%
% op_roty - returns the matrix corresponding to a counter clockwise rotation
% of "angle" degrees around the y axis.
%
% Required input variables:
% angle: rotation angle in degrees.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@fluor.quimica.uniovi.es>
% Created: June 2011

c = cos(angle*pi/180);
s = sin(angle*pi/180);
m = [c 0 s; 0 1 0; -s 0 c];

endfunction
