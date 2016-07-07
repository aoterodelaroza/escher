% Copyright (c) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function g = grid_readcube (cubefile, LOG=0)
% function g = grid_readcube (cube, LOG=0)
%
% grid_readcube - read a grid from a gaussian cube file.
%
% Required input variables:
% {cubefile}: name of the cube file.
%
% Output:
% g: grid object.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: Jan 2012

  bohrtoans = 0.52917720859;

  [fqub,msg] = fopen(cubefile,'r');
  if (fqub < 0 || ferror(fqub))
    disp(msg)
    error("mol_readcube: Could not open -- %s",cubefile);
  endif

  g = grid_();

  ## title lines in the cube file
  g.name = fgetl(fqub); title{2} = fgetl(fqub);

  ## read in number of atoms, origin, and grid
  line = fgetl(fqub);
  [tok,line] = strtok(line); nat = str2num(tok);
  g.x0 = str2num(line);
  g.a = g.dx = zeros(3,3);
  for i = 1:3
    line = fgetl(fqub);
    [tok,line] = strtok(line); g.n(i) = str2num(tok);
    g.dx(i,:) = str2num(line);
    g.a(i,:) = g.dx(i,:) * g.n(i);
  endfor
  g.omega = det(g.a);

  ## read in atoms and atomic positions
  for i = 1:nat
    line = fgetl(fqub);
  endfor

  ## now read in the function values at the grid
  alldata = fscanf(fqub,"%g");
  fclose(fqub);

  ## rearrange the grid
  nd = 0;
  for ix = 1:g.n(1)
    for iy = 1:g.n(2)
      g.f(ix,iy,1:g.n(3)) = alldata(nd+1:nd+g.n(3));
      nd = nd + g.n(3);
    endfor
  endfor

  if (LOG>0)
    printf('grid_readcube:\n');
    printf('File read: %s\n', cubefile);
    printf('Grid origin: %.6f %.6f %.6f\n', g.x0);
    printf('Grid dimensions: %d %d %d\n', g.n);
    printf('Grid matrix:\n');
    for i = 1:3
      printf('%15.9f %15.9f %15.9f\n', g.dx(i,1:3));
    endfor
    printf('Number of data values read: %d\n', nd);
    printf('Shape of cubedata matrix: %d %d %d\n', size(g.f));
    printf('int(f): %.7f \n', sum(sum(sum(g.f))) * g.omega / prod(g.n));
    printf('int(f^2): %.7f \n', sum(sum(sum(g.f.^2))) * g.omega / prod(g.n));
    printf('int(|f|): %.7f \n', sum(sum(sum(abs(g.f)))) * g.omega / prod(g.n));
  endif

endfunction
