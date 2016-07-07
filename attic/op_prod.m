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

function T = op_prod (R, S)
% function T = op_prod (R, S)
%
% op_prod - product of two rotation-translation operations (Seitz op.).
% Each operation is 3x4 and it is assumed to work on a 3x1 column vector:
%
%       R S x = R (S x) = (RS) x = T x
%
% Required input variables:
% R,S: 3x4 rotation-traslation matrices.
%
% Optional input variables (all with default values):
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: July 2011

T(1:3,1:3) = R(1:3,1:3) * S(1:3,1:3);
T(:,4) = R(1:3,1:3) * S(1:3,4) + R(1:3,4);

endfunction
