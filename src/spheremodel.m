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

function [vv ifv] = spheremodel(ndiv=1)
% function [vv ifv] = spheremodel(ndiv=1)
%
% spheremodel - create a vertex and faces polyhedral representation 
%               of a sphere. Center is (0,0,0) and all vectors normalized.
%               
%
% Optional input variables:
% ndiv: number of subdivision of the initial icosahedron.
% 
% Output variables:
% vv(nv,3): list of vertices.
% ifv(nf,3): integer list of faces (first vertex is 1).
%

  ## from tessel: icosahedron as 3 intersecting golden-ratio rectangles. 
  tau = (1 + sqrt(5))/2;
  phi = tau - 1;
  rad0 = sqrt(3 - tau);

  vv = [
        1, 0,  phi
        1, 0, -phi
       -1, 0,  phi
       -1, 0, -phi
        0,  phi ,1
        0, -phi ,1
        0,  phi ,-1
        0, -phi ,-1
         phi,  1, 0
        -phi,  1, 0
         phi, -1, 0
        -phi, -1, 0
        ] / rad0;

  ## faces
  ifv = [
         1,     9,    5
         10,    9,    7
         3,     5,   10
         9,    10,    5
         3,   10,    4
         10,    7,    4
         1,   11,    2
         9,    1,    2
         7,    9,    2
         12,    3,    4
         6,    1,    5
         11,    1,    6
         12,   11,    6
         5,    3,    6
         6,    3,    12
         8,    7,    2
         8,    2,   11
         8,   11,   12
         4,    7,    8
         12,    4,    8
         ];
  nvo = nv = size(vv,1); nf = size(ifv,1);

  for i = 1:ndiv
    ## allocate space for the new vertices
    vv(nv+1:nv+nf*3/2,1:3) = 0;
    ifv(nf+1:4*nf,1:3) = 0;

    ## create new vertices and faces
    icon = zeros(nv);
    for j = 1:nf
      ia = ifv(j,1); ib = ifv(j,2);
      if (icon(ia,ib) == 0)
        nv = nv + 1;
        vv(nv,:) = (vv(ia,:) + vv(ib,:)) / 2;
        icon(ib,ia) = icon(ia,ib) = in1 = nv;
      else
        in1 = icon(ia,ib);
      endif
      ia = ifv(j,2); ib = ifv(j,3);
      if (icon(ia,ib) == 0)
        nv = nv + 1;
        vv(nv,:) = (vv(ia,:) + vv(ib,:)) / 2;
        icon(ib,ia) = icon(ia,ib) = in2 = nv;
      else
        in2 = icon(ia,ib);
      endif
      ia = ifv(j,3); ib = ifv(j,1);
      if (icon(ia,ib) == 0)
        nv = nv + 1;
        vv(nv,:) = (vv(ia,:) + vv(ib,:)) / 2;
        icon(ib,ia) = icon(ia,ib) = in3 = nv;
      else
        in3 = icon(ia,ib);
      endif
      of = ifv(j,:);
      ifv(j,:) = [of(1) in1 in3];
      nf = nf + 1; ifv(nf,:) = [of(2) in2 in1];
      nf = nf + 1; ifv(nf,:) = [of(3) in3 in2];
      nf = nf + 1; ifv(nf,:) = [in1 in2 in3];
    endfor

    ## renormalize
    for i = nvo+1:nv
      vv(i,:) /= norm(vv(i,:));
    endfor
  endfor

endfunction
