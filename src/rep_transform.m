% Copyright (C) 2011--12 Victor Lua~na and Alberto Otero-de-la-Roza
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

function rep = rep_transform(repi, rot=[1 0 0; 0 1 0; 0 0 1], tr=[0 0 0]);
% function rep = rep_transform(repi, rot=[1 0 0; 0 1 0; 0 0 1], tr=[0; 0; 0]);
%
% rep_transform - rotate and translate a representation
%
% Required input variables:
% {repi}: input representation
%
% Optional input variables:
% rot: rotation matrix (3x3). x_trans = x_orig * rot + tr, with x row vectors.
% tr: translation vector (3)
%
% Required output variables:
% {rep}: output representation.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

  rep = repi;
  if (all(size(tr) == [3 1]))
    tr = tr';
  endif

  for i = 1:rep.nball
    rep.ball{i}.x = repi.ball{i}.x * rot + tr;
  endfor
  for i = 1:rep.nstick
    rep.stick{i}.x0 = repi.stick{i}.x0 * rot + tr;
    rep.stick{i}.x1 = repi.stick{i}.x1 * rot + tr;
  endfor
  for i = 1:rep.nvertex
    rep.vertex{i}.x = repi.vertex{i}.x * rot + tr;
  endfor

endfunction
