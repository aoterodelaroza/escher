% Copyright (C) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function err = mol_writeg09(mol, code="def", conditions=struct([]), LOG=2)
%function err = mol_writeg09(mol, code="def", conditions, LOG=0)
%
% mol_writeg09 - write a gaussian09 input file.
%
% Required input variables:
% mol: description of the molecule.
%
% Optional input variables (all have default values):
% {code="def"): Root for the filenames (.chk, .com, .wfn, ...).
%                     By defult, the standard output is used.
% {conditions}: modification of the default calculation conditions.
%     The default options are:
%     *.method = "RHF",
%     *.basis = "3/21g",
%     *.opt = "",
%     *.doit = 0  (change to 1 to effectively do the calculation)
% {calc}: the "method/basis_set" description of the calculation.
% {opt}: single point, optimization, frequency, ... The gaussian task to
%        perform.
% {LOG = 1}: print information about the data read in if LOG>0.
%
% Required output variables:
% err: error code. 0 means no error.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: January 2012
% ReCreated: January 2013

  angtobohr = 1.88972613288564;
  bohrtoans = 0.52917720859;
  err = 0;

% Default Conditions:

  cdn = struct(
      "basis", "aug-cc-PVTZ",
      "method", "b3lyp",
      "opt", "opt=(tight,Newton,CalcAll)",
      "freq", [],
      "doit", 0
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

  if (isfield(mol,"name"))
    title = mol.name;
  else
    title = "molecule";
  endif

  err = 1;
  filein  = sprintf("%s.com", code);
  fileout = sprintf("%s.log", code);
  filewfn = sprintf("%s.wfn", code);
  filechk = sprintf("%s.chk", code);
  [fid,msg] = fopen(filein,"w+");
  if (fid < 0 || ferror(fid)) 
    disp(msg)
    printf("mol_writeg09: Could not find or create -- %s",filein);
    err=1;
    return
  endif

  calc = sprintf("%s/%s %s %s", cdn.method, cdn.basis, cdn.opt, cdn.freq);
  fprintf(fid, "$RunGauss\n");
  fprintf(fid, "%%mem=2GB\n");
  fprintf(fid, "%%NoSave\n");
  fprintf(fid, "%%NProcShared=4\n");
  fprintf(fid, "%%chk=%s\n", filechk);
  fprintf(fid, "#P %s %s Int=(Grid=UltraFine) scf=tight 5d 7f\n", calc);
  fprintf(fid, "gfinput gfoldprint pop=(full,esp) density=current output=wfn\n");
  fprintf(fid, "IOP(6/7=3)\n");
  fprintf(fid, "\n");
  fprintf(fid, "(VLC00) %s (%s)\n", title, calc);
  fprintf(fid, "\n");
  fprintf(fid, "0 1\n");
  for i = 1 : mol.nat;
     fprintf(fid,"   %-2s  ", mol.atname{i});
     fprintf(fid,"  %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
     fprintf(fid,"\n");
  endfor
  fprintf(fid, "\n");
  fprintf(fid, "%s\n\n", filewfn);

  fclose(fid);

  if (LOG > 0)
     printf("mol_writeg09: g09 calculation\n");
     printf("Created file input: %s\n", filein);
     for [val, key] = conditions
        printf("%s: %s\n", key, val);
     endfor
     if (LOG > 1)
        printf("Molecule: %s\n", mol.name);
        for i = 1 : mol.nat;
           printf("%5d %-2s %5d", mol.atname{i}, mol.atnumber(i));
           printf("  %15.9f %15.9f %15.9f", mol.atxyz(1:3,i));
           printf("\n");
        endfor
     endif
  endif

  if (cdn.doit == 1)
     order = sprintf("g09 < %s > %s", filein, fileout);
     system(order);
  endif

endfunction
