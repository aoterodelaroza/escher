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

function [mask] = mask_merge(varargin);
% function [mask] = mask_merge(mask1, mask2,...)
%
% mask_merge - merge masks.
%
% Required input variables:
% {mask1,mask2,...}: masks or cell arrays of masks to merge.
%
% Required output variables:
% {mask}: merged mask. 
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: December 2011

mask = crmask();
n = 0;
for i = 1:length(varargin)
  if (!iscell(varargin{i}))
    mcell = {varargin{i}};
  else
    mcell = varargin{i};
  endif

  for k = 1:length(mcell)
    mask1 = mcell{k};
    n1 = mask1.nat;
    for j = 1 : n1
      n++;
      mask.l(n,1:3) = mask1.l(j,1:3);
      mask.i(n) = mask1.i(j);
    endfor
  endfor
endfor
mask.nat = n;

endfunction
