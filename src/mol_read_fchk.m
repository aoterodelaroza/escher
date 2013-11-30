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

function [mol] = mol_read_fchk(rootfile, LOG=1)
%function [mol] = mol_read_fchk(rootfile, LOG=1)
%
% mol_read_fchk - read in a mol description from a gaussian g09
% fchk file.
%
% Required input variables:
% {rootfile}: Add ".fchk" and this will be the input file.
% Atom coordinates are assumed to be in bohr.
%
% Optional input variables (all have default values):
% {LOG}: print the final result if LOG>0.
%
% Ouput variables:
% {mol}: struct with the mol description.
%

   bohrtoans = 0.52917720859;

   fchkfile = sprintf("%s.fchk", rootfile);
   if (!exist(fchkfile))
      error(sprintf("Could not find file: %s\n",fchkfile));
   endif
   fid = fopen(fchkfile,"r");

   mol=molecule();

   nline=0;
   again = 1;
   while (again) 
      ++nline;
      line = fgetl(fid);
      if(nline==1)
         mol.name = line;
      elseif(strfind(line,"Number of atoms"))
         [d1 d2 d3 d4 mol.nat] = sscanf(line, "%s %s %s %s %d","C");
      elseif(strfind(line,"Atomic numbers"))
         line = fgetl(fid);
         iat=0;
         while (iat<mol.nat && !feof(fid))
            [li, nl] = sscanf(line, "%d", Inf);
            for i = 1:nl
               if (i <= mol.nat)
                  mol.atnumber(++iat) = li(i);
                  [atsymb,atprop] = mol_dbsymbol(mol.atnumber(iat));
                  mol.atmass(iat) = atprop.mass;
                  mol.atname{iat} = atsymb;
               endif
            endfor
            if (iat < mol.nat)
               line = fgetl(fid);
            endif
         endwhile
      elseif(strfind(line,"Current cartesian coordinates"))
         mol.atxyz=zeros(3,mol.nat);
         x = zeros(1,3*mol.nat);
         line = fgetl(fid);
         n = 3*mol.nat;
         i = 0;
         while (i<n && !feof(fid))
            [xx,nnx] = sscanf(line, "%f", Inf);
            for j = 1 : nnx
               x(++i) = xx(j);
            endfor
            if (i < n)
               line = fgetl(fid);
            endif
         endwhile
         mol.atxyz = reshape(x,3,mol.nat);
         %k = 0;
         %for i = 1:mol.nat
         %   for j = 1:3
         %      mol.atxyz(j,i) = x(++k);
         %   endfor
         %endfor
         mol.atxyz = mol.atxyz * bohrtoans;
      elseif(strfind(line,"Force Field"))
         again = 0;
      endif
   endwhile

   if (LOG > 0)
      printf("File.: %s\n", fchkfile);
      printf("Title: %s\n", mol.name);
      printf("--i- -zi- symb --------Coord (Angfstrom)-------\n");
      for i = 1:mol.nat
         zi = mol.atnumber(i);
         symb = mol_dbsymbol(zi);
         printf("%4d %4d %2s", i, zi, symb);
         printf("%12.6f %12.6f %12.6f\n", mol.atxyz(:,i));
      endfor
   endif

endfunction
