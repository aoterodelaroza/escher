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

function mplot_writepov(rootname, runpov=1, size=[600,600], LOG=1))
% function mplot_writepov(rootname="mol", runpov=1, size=[600,600], LOG=1))
%
% mplot_writepov - create the input for povray and, eventually, run the
%    povray and convert programs to produce a PNG file.
%
% Required input variables: NONE
%
% Optional input variables (all have default values):
% {rootname}: Use this as root for the ".pov", ".tga", and ".png" files that
%   can be created by this routine.
% {runpov}: povray and convert will be run only if this switch is set to true.
% {size}: size of the final image.
% {LOG}: print the final result if LOG>0.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

global ge 
global pov

filepov = sprint("%s.%s", rootname, '.pov');
filetga = sprint("%s.%s", rootname, '.tga');
filepng = sprint("%s.%s", rootname, '.png');

[fpov,msg] = fopen(filepov,'w+');
if (fpov < 0 || ferror(fpov))
  disp(msg)
  error("mplot_writepov: Could not open -- %s",filepov);
endif

# ge variables:
ge_maxtrace = 10;

# Get limits for points in the graphical database:
xmin = min(ge.pt')'; xmax = max(ge.pt')';
xct = (xmax + xmin)/2; xdel = xmax - xmin;
printf(fpov,'#include "colors.inc"\n');
printf(fpov,'#include "textures.inc"\n');
printf(fpov,'#include "woods.inc"\n');
printf(fpov,'#include "stones.inc"\n');
printf(fpov,'// Model created by MolWare\n');
printf(fpov,'// (c) 2011 Victor Lua~na and Alberto Otero-de-la-Roza\n');
printf(fpov,'// Departamento de Quimica Fisica, Universidad de Oviedo\n');
printf(fpov,'// 33006-Oviedo, Spain\n');
printf(fpov,'// victor@carbono.quimica.uniovi.es\n');
printf(fpov,'#max_trace_level %d\n', ge_maxtrace);

fclose(fpov);


endfunction
