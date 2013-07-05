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

function [vv ifv] = cylindermodel(nbase=16)
% function [vv ifv] = cylindermodel(ndiv=0)
%
% cylindermodel - create a vertex and faces polyhedral representation 
%                 of a cylinder. The center of the bases are at (0,0,0)
%                 and (0,0,1) and their radius is 1.
%               
%
% Optional input variables:
% nbase: number of vertices in each base.
% 
% Output variables:
% vv(nv,3): list of vertices.
% ifv(nf,3): integer list of faces (first vertex is 1).
%

  ## polygon vertices
  vv = zeros(2*nbase+2,3);
  z = i; sz = cos(2*pi/nbase) + i * sin(2*pi/nbase);
  for i = 1:nbase
    z = z * sz;
    vv(i,:) = [real(z) imag(z) 0];
    vv(nbase+i,:) = [real(z) imag(z) 1];
  endfor
  vv(2*nbase+1,:) = [0 0 0];
  vv(2*nbase+2,:) = [0 0 1];

  ## polygon faces
  ifv = zeros(4*nbase,3);
  nf = 0;
  for i = 1:nbase
    ic1 = i; ic2 = mod(i,nbase)+1;
    nf = nf + 1; ifv(nf,:) = [2*nbase+1, ic1, ic2]; # triangle on first base
    nf = nf + 1; ifv(nf,:) = [2*nbase+2, nbase+ic1, nbase+ic2]; # triangle on second base
    nf = nf + 1; ifv(nf,:) = [ic1, ic2, nbase+ic1]; # tube triangle 1
    nf = nf + 1; ifv(nf,:) = [nbase+ic2, nbase+ic1, ic2]; # tube triangle 2
  endfor

endfunction
