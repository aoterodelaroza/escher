#! /usr/bin/octave -q
% Copyright (C) 2013 Victor Lua~na and Alberto Otero-de-la-Roza
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

function err = mol_fchk2topo (name, conditions=struct([]), LOG=1)
% function err = mol_fchk2topo (name, conditions=struct([]), LOG=1)
%
% mol_fchk2topo - read a g09 fchk file and determine the MEP/rho topology
% using cubegen and critic2g.
%
% The predetermined sequence of calculations is:
% 1. Read the fchk file and determine the range for cubegen
% 2. Get the grid description of rho and MEP using cubegen
% 3. Run critic2g to analyze the topologies of both fields
%
% Required input variables:
% {name}: rootname for the next files.
%    <name>.fchk ------> fchk file used by cubegen and this routine.
%
% Files generated by the run of cubegen and critic2g:
%    <name>-rho.cube --> cube files for rho and mep.
%    <name>-mep.cube
%    <name>-rho.crin --> input to critic2g.
%    <name>-rho.crout
%    <name>-mep.crout
%    <name>-mep.crout
%
% Optional input variables (all have default values):
% {conditions}: modification of the default calculation conditions.
%    The default options are defined below, under the cdn struct:
%    *.qal (6) quality level of the rho/mep grids.
%    *.mrg (2) empty margin surrounding the finite molecule.
%    *.rho: topology of rho? (0/1)
%    *.mep: topology of mep? (0/1)
%    *.cube: a cube rather than a parallepipedic grid? (0/1)
%    *.cytpe}: type of gaussian calculation. Options: SCF, MP2, CI, QCI.
% {LOG=1}: print information about the data read in if LOG>0.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: Jan 2013

   bohrtoans = 0.52917720859;

% Default conditions:
   cdn = struct(
      "qal", 6,
      "mrg", 2,
      "rho", 1,
      "mep", 1,
      "cube", 0,
      "ctype", "scf"
      );

% Should the default conditions change?

   for [val, key] = conditions
      if (isfield(cdn,key))
         cdn.(key) = val;
         # cdn = setfield (key, val);
      else
         printf("mol_writeg09: Unknown option %s = %s", key, val);
         err = 1;
         error("error found!");
         return
      endif
   endfor

% Time the runtime
   tic;

% Read the fchk y determine the size of the cube grid
   mol = mol_read_fchk(name,0);
   mol.atxyz = mol.atxyz * bohrtoans;
   [x0,nx,dx] = mol_cuberange(mol, cdn.qal, cdn.mrg, cdn.cube, 0);

% Create a bash script to run cubegen and critic2g
% (this is required by the form cubegen is used)
   fsh = sprintf("%s-cr2.sh", name);
   [fid,msg] = fopen(fsh, "w+");
   if (fid < 0 || ferror(fid))
      disp(msg)
      error("mol_fchk2topo: could not open -- %s", fsh);
   endif

   if (cdn.rho > 0)
      fprintf(fid, "cat << EOF | ");
      fprintf(fid, "cubegen 0 density=%s %s.fchk %s-rho.cube -1 h\n", \
              cdn.ctype, name, name);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", -6, x0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(1), dx(1), 0.0, 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(2), 0.0, dx(2), 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(3), 0.0, 0.0, dx(3));
      fprintf(fid, "EOF\n\n");
      fprintf(fid, "critic2g << EOF > %s-rho.crout\n", name);
      fprintf(fid, "nohmove\n");
      fprintf(fid, "crystal\n");
      fprintf(fid, "   spg P 1\n");
      fprintf(fid, "   struct %s-rho.cube cube\n", name);
      fprintf(fid, "   den %s-rho.cube cube\n", name);
      fprintf(fid, "   nocore\n", name);
      fprintf(fid, "endcrystal\n");
      fprintf(fid, "auto\n");
      fprintf(fid, "EOF\n");
      fprintf(fid, "mv stdin_cps.xyz %s_rho.xyz\n\n",name);
   endif

   if (cdn.mep > 0)
      fprintf(fid, "cat << EOF | ");
      fprintf(fid, "cubegen 0 potential=%s %s.fchk %s-mep.cube -1 h\n", \
              cdn.ctype, name, name);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", -6, x0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(1), dx(1), 0.0, 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(2), 0.0, dx(2), 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(3), 0.0, 0.0, dx(3));
      fprintf(fid, "EOF\n\n");
      fprintf(fid, "critic2g << EOF > %s-mep.crout\n", name);
      fprintf(fid, "crystal\n");
      fprintf(fid, "   spg P 1\n");
      fprintf(fid, "   struct %s-mep.cube cube\n", name);
      fprintf(fid, "   den %s-mep.cube cube\n", name);
      fprintf(fid, "   nocore\n", name);
      fprintf(fid, "endcrystal\n");
      fprintf(fid, "auto\n");
      fprintf(fid, "EOF\n");
      fprintf(fid, "mv stdin_cps.xyz %s_mep.xyz\n\n",name);
   endif

   fclose(fid);

   if (cdn.rho > 0 || cdn.mep > 0)
      system(sprintf("sh %s", fsh));
   endif

   cpu = toc;

   if (LOG > 0)
      printf("Grid:\n");
      printf("%5d %12.6f %12.6f %12.6f\n", -6, x0);
      printf("%5d %12.6f %12.6f %12.6f\n", nx(1), dx(1), 0.0, 0.0);
      printf("%5d %12.6f %12.6f %12.6f\n", nx(2), 0.0, dx(2), 0.0);
      printf("%5d %12.6f %12.6f %12.6f\n", nx(3), 0.0, 0.0, dx(3));
      if (cdn.rho > 0)
         printf("Check topo(rho) in %s-rho.crout  [%s]\n", name, cdn.ctype); 
      endif
      if (cdn.mep > 0)
         printf("Check topo(rho) in %s-mep.crout  [%s]\n", name, cdn.ctype); 
      endif
      printf("CPU time: %.4f\n", cpu);
   endif

endfunction

# Use it as a script
if (!exist("argn"))
  if (nargin > 0)
    args = argv();
    err = mol_fchk2topo(args{1});
  else
    error("Use it as mol_fchk2topo <fchk_file> (skip the .fchk)!")
  endif
endif
