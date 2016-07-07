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

function [mol] = mol_smiles2topo (smiles="O=C(O)C(N)C Ala", code="mol", conditions=[], LOG=1)
% function mol_smiles2topo (smiles, code, conditions, LOG=1)
%
% mol_smiles2topo - Starting from the smiles code or a molecule, this rutine
% performs a complex collection of calculations that produce the topological
% analysis of the molecular electrostatic potential. This routine has been
% initially designed to determine the lone pairs of the aminoacids.
%
% Required input variables:
% {smiles}:
% {code}: root for the project filenames. Mainly:
%    <code>.xyz ------> starting cartesian coordinates.
%    <code>.fchk -----> fchk file used by cubegen and this routine.
%    <code>-rho.cube -> cube files for rho and mep.
%    <code>-mep.cube
%    <code>-cr2.sh ---> sh script to run the calculations
% {conditions}: structure containing variations of the default conditions.
%
% Optional input variables (all have default values):
% {LOG = 1}: print information about the data read in if LOG>0.
%            LOG = 0  no output.
%            LOG = 1  debug information.
%
% Final files containing relevant data:
% 
% 
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: Jan 2013

   tic;

   bohrtoans = 0.52917720859;

% Default Conditions:

   cdn = struct(
      "basis", "aug-cc-PVTZ",
      "method", "b3lyp",
      "opt", "opt=(tight,Newton,CalcAll)",
      "rho", 0,
      "mep", 1,
      "qual", 6,
      "marg", 3,
      "cube", 1,
      "ctype", "scf"
      );

% Should the default conditions change?

   for [val, key] = conditions
      cdn = setfield (key, val);
      if (isfield(orig,key))
         orig.(key) = val;
      endif
   endfor

% Create and read an xyz struture from the smiles:

   fxyz = sprintf("%s.xyz", code);
   order = sprintf("obabel -:'%s' -oxyz -h --gen3d -O %s", smiles, fxyz);
   system(order);

   mol = mol_readxyz(fxyz);

% Prepare and perform the g09 calculation:

   calc = sprintf("%s/%s", cdn.method, cdn.basis);
   opt = sprintf("%s", cdn.opt);
   err = mol_writeg09(mol, code, calc, opt);
   system(sprintf("formchk %s.chk", code));

   return

% Prepare and perform the cubegen calculation:

   mol = mol_readfchk(code,0);
   mol.atxyz = mol.atxyz * bohrtoans;
   [x0,nx,dx] = mol_cuberange(mol, cdn.qal, cdn.mrg, cdn.cube, 0);

   fsh = sprintf("%s-cr2.sh", code);
   [fid,msg] = fopen(fsh, "w+");
   if (fid < 0 || ferror(fid))
      disp(msg)
      error("mol_cr2topo: could not open -- %s", fsh);
   endif

   if (cdn.rho > 0)
      fprintf(fid, "cat << EOF | ");
      fprintf(fid, "cubegen 0 density=%s %s.fchk %s-rho.cube -1 h\n", \
              ctype, code, code);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", -6, x0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(1), dx(1), 0.0, 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(2), 0.0, dx(2), 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(3), 0.0, 0.0, dx(3));
      fprintf(fid, "EOF\n\n");
      fprintf(fid, "critic2g << EOF > %s-rho.cr2out\n", code);
      fprintf(fid, "nohmove\n");
      fprintf(fid, "crystal\n");
      fprintf(fid, "   spg P 1\n");
      fprintf(fid, "   struct %s-rho.cube cube\n", code);
      fprintf(fid, "   den %s-rho.cube cube\n", code);
      fprintf(fid, "   nocore\n", code);
      fprintf(fid, "endcrystal\n");
      fprintf(fid, "auto\n");
      fprintf(fid, "EOF\n");
      fprintf(fid, "mv stdin_cps.xyz %s_rho.xyz\n\n",code);
   endif

   if (cdn.mep > 0)
      fprintf(fid, "cat << EOF | ");
      fprintf(fid, "cubegen 0 potential=%s %s.fchk %s-mep.cube -1 h\n", \
              ctype, code, code);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", -6, x0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(1), dx(1), 0.0, 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(2), 0.0, dx(2), 0.0);
      fprintf(fid, "%5d %12.6f %12.6f %12.6f\n", nx(3), 0.0, 0.0, dx(3));
      fprintf(fid, "EOF\n\n");
      fprintf(fid, "critic2g << EOF > %s-mep.cr2out\n", code);
      fprintf(fid, "crystal\n");
      fprintf(fid, "   spg P 1\n");
      fprintf(fid, "   struct %s-mep.cube cube\n", code);
      fprintf(fid, "   den %s-mep.cube cube\n", code);
      fprintf(fid, "   nocore\n", code);
      fprintf(fid, "endcrystal\n");
      fprintf(fid, "auto\n");
      fprintf(fid, "EOF\n");
      fprintf(fid, "mv stdin_cps.xyz %s_mep.xyz\n\n",code);
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
      if (rho > 0)
         printf("Check topo(rho) in %s-rho.cr2out  [%s]\n", code, ctype); 
      endif
      if (mep > 0)
         printf("Check topo(rho) in %s-mep.cr2out  [%s]\n", code, ctype); 
      endif
      printf("CPU time: %.4f\n", cpu);
   endif

endfunction
