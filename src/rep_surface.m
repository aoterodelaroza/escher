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

function rep = rep_surface(addto="",f,uv0=[-1 -1],uv1=[1 1],nuv=[41 41],closed=[0 0],...
                           rgb0=[255 0 0],tex="opaque_triangle_default",...
                           grid=0,grrad=0.005,grrgb=[255 0 0],grtex="stick_default")
% function rep = rep_surface(addto="",f,uv0=[-1 -1],uv1=[1 1],nuv=[41 41],closed=[0 0],...
%                            rgb0=[255 0 0],tex="opaque_triangle_default",...
%                            grid="",grrad=0.005,grrgb=[255 0 0],grtex="stick_default")
%
% rep_surface - surface
%
% Required input variables:
%
% Optional input variables (all have default values):
%
  
  ## uv-step
  ud = (uv1 - uv0) ./ nuv;
  isup = nuv-1;
  if (closed(1) != 0)
    ud(1) = (uv1(1) - uv0(1)) / (nuv(1)-1);
    isup(1) = nuv(1);
  endif
  if (closed(2) != 0)
    ud(2) = (uv1(2) - uv0(2)) / (nuv(2)-1);
    isup(2) = nuv(2);
  endif

  ## initial representation 
  if (!isempty(addto) && isstruct(addto))
    rep = addto;
  else
    rep = representation();
  endif

  ## calculate function values
  u = zeros(nuv(1)*nuv(2),2);
  n = 0;
  for i = 1:nuv(1)
    for j = 1:nuv(2)
      n++;
      u(n,:) = uv0 + ([i j]-1) .* ud;
    endfor
  endfor
  fval = f(u);

  ## calcualte color values
  if (isnumeric(rgb0))
    rgb = zeros(nuv(1)*nuv(2),5);
    rgb0 = fillrgb(rgb0);
    for i = 1:nuv(1)*nuv(2)
      rgb(i,:) = rgb0;
    endfor
  elseif (ischar(rgb0))
    rgb = zeros(nuv(1)*nuv(2),5);
    rgb0 = fillrgb(color(rgb0));
    for i = 1:nuv(1)*nuv(2)
      rgb(i,:) = rgb0;
    endfor
  else
    rgb = fillrgb(rgb0(u));
  endif

  ## fill coordinates and colors
  n = 0;
  nv0 = rep.nvertex;
  for i = 1:nuv(1)
    for j = 1:nuv(2)
      ## Coordinates
      n++;
      rep.vertex{n} = vertex();
      rep.vertex{n}.x = fval(n,:);
      rep.vertex{n}.rgb = rgb(n,:);
    endfor
  endfor
  rep.nvertex += n;

  ## add triangles
  [rep itex] = rep_registertexture(rep,tex);
  if (grid > 0)
    [rep igrtex] = rep_registertexture(rep,grtex);
  endif
  n = rep.ntriangle;
  m = rep.nstick;
  for i = 1:isup(1)
    i1 = mod(i,nuv(1)) + 1;
    for j = 1:isup(2)
      j1 = mod(j,nuv(2)) + 1;

      idx = [(i-1)*nuv(2)+j, (i1-1)*nuv(2)+j, (i1-1)*nuv(2)+j1, (i-1)*nuv(2)+j1];

      n++;
      rep.triangle{n} = triangle();
      rep.triangle{n}.idx = nv0 + [idx(1) idx(2) idx(4)];
      rep.triangle{n}.rgb = (rep.vertex{rep.triangle{n}.idx(1)}.rgb + ...
                             rep.vertex{rep.triangle{n}.idx(2)}.rgb + ...
                             rep.vertex{rep.triangle{n}.idx(3)}.rgb)/3;
      rep.triangle{n}.tex = itex;
      n++;
      rep.triangle{n} = triangle();
      rep.triangle{n}.idx = nv0 + [idx(4) idx(2) idx(3)];
      rep.triangle{n}.rgb = (rep.vertex{rep.triangle{n}.idx(1)}.rgb + ...
                             rep.vertex{rep.triangle{n}.idx(2)}.rgb + ...
                             rep.vertex{rep.triangle{n}.idx(3)}.rgb)/3;
      rep.triangle{n}.tex = itex;

      if (grid > 0)
        m++;
        rep.stick{m} = stick();
        rep.stick{m}.name = "surface";
        rep.stick{m}.x0 = fval(idx(1),:);
        rep.stick{m}.x1 = fval(idx(2),:);
        rep.stick{m}.r = grrad;
        rep.stick{m}.rgb = grrgb;
        rep.stick{m}.tex = igrtex;
        m++;
        rep.stick{m} = stick();
        rep.stick{m}.name = "surface";
        rep.stick{m}.x0 = fval(idx(1),:);
        rep.stick{m}.x1 = fval(idx(4),:);
        rep.stick{m}.r = grrad;
        rep.stick{m}.rgb = grrgb;
        rep.stick{m}.tex = igrtex;
        m++;
        rep.stick{m} = stick();
        rep.stick{m}.name = "surface";
        rep.stick{m}.x0 = fval(idx(2),:);
        rep.stick{m}.x1 = fval(idx(4),:);
        rep.stick{m}.r = grrad;
        rep.stick{m}.rgb = grrgb;
        rep.stick{m}.tex = igrtex;
      endif
    endfor
  endfor
  rep.ntriangle = n;
  rep.nstick = m;

endfunction

