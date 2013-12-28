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

function [rep] = rep_merge(varargin);
% function [rep] = rep_merge(rep1, rep2,...)
%
% rep_merge - merge representations. The camera and background color are
% taken from the last representation that has some.
%
% Required input variables:
% {rep1,rep2,...}: representations or cell arrays of representations to join.
%
% Required output variables:
% {rep}: merged representation.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

  rep = representation();
  rep.name = "merged: ";

  n = 0;
  for i = 1:length(varargin)
    if (!iscell(varargin{i}))
      mcell = {varargin{i}};
    else
      mcell = varargin{i};
    endif

    for k = 1:length(mcell)
      rep1 = mcell{k};
      if (isfield(rep1,"name"))
        rep.name = sprintf("%s %s",rep.name,rep1.name);
      endif

      ## balls
      n1 = rep1.nball; 
      for j = 1:n1
        rep.nball += 1;
        rep.ball{rep.nball} = rep1.ball{j};
      endfor

      ## sticks
      n1 = rep1.nstick; 
      for j = 1:n1
        rep.nstick += 1;
        rep.stick{rep.nstick} = rep1.stick{j};
      endfor

      ## vertices
      nv0 = rep.nvertex;
      n1 = rep1.nvertex; 
      for j = 1:n1
        rep.nvertex += 1;
        rep.vertex{rep.nvertex} = rep1.vertex{j};
      endfor

      ## triangles
      n1 = rep1.ntriangle; 
      for j = 1:n1
        rep.ntriangle += 1;
        rep.triangle{rep.ntriangle} = rep1.triangle{j};
        rep.triangle{rep.ntriangle}.idx += nv0;
      endfor

      ## lights
      n1 = rep1.nlight; 
      for j = 1:n1
        rep.nlight += 1;
        rep.light{rep.nlight} = rep1.light{j};
      endfor

      ## lights
      n1 = rep1.nsurf; 
      for j = 1:n1
        rep.nsurf += 1;
        rep.surf{rep.nsurf} = rep1.surf{j};
      endfor

      ## camera and background color
      rep.cam = rep1.cam;
      rep.bgcolor = rep1.bgcolor;
    endfor
  endfor

endfunction
