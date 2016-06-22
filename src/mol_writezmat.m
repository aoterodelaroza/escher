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

function xzmat = mol_writezmat(mol,file="",header=0,izmat=[],izvar=[])
% function xzmat = mol_writezmat(mol,file="",header=0,izmat=[],izvar=[])
%
% mol_writezmat - write a zmatrix calculated from the Cartesian
%   coordinates in mol. If no izmat is present, the atoms are taken in the
%   same sequence as in mol. If izmat is present, it should have the structure:
%   izmat = [ 
%             i1 i2 i3 i4 
%             i5 i6 i7 i8
%             ...
%           ]
%   where atom i1 in mol should appear in the first z-matrix entry,
%   and it connects to i2 (distance), i3 (angle), and i4
%   (dihedral). The 2:4 entries in the first line, 3:4 in the second line, and 
%   fourth entry in the third line are ignored.
%
%   If izvar is present, make some of the entries in the z-matrix "variables"
%   in Gaussian's input format. izvar should have the same structure as izmat.
%   If an entry in izvar is non-zero, it is made a variable.
% 
% Required input variables:
% mol: the input molecule
%
% Optional input variables:
% file: output file (zmat). If not present, inhibit output.
% header: if true (1), write a proper Gaussian input file.
% izmat: array of atomic indices for the z-matrix construction. 
% If not present, use the atomic sequence in input order.
%
% Output variables:
% xzmat: array containing distances, angles, and dihedrals.
%

  ## open the output file
  if (!isempty(file))
    fid = fopen(file,"w");
  else
    fid = fopen("/dev/null","w");
  endif

  ## fill the izmat
  if (isempty(izmat)) 
    izmat = zeros(mol.nat,4);
    for i = 1:mol.nat
      izmat(i,:) = max([i i-1 i-2 i-3],0);
    endfor
  else
    izmat(1,2:4) = 0;
    izmat(2,3:4) = 0;
    izmat(3,4) = 0;
  endif
  
  ## write the header, if applicable
  if (header)
    nelec = sum(mol.atnumber);
    if (mod(nelec,2) == 0)
      mult = 1;
    else
      mult = 2;
    endif
    fprintf(fid, "%%mem=2GB\n");
    fprintf(fid, "%%nproc=4\n");
    fprintf(fid, "#p b3lyp sto-3g\n");
    fprintf(fid, "\n");
    fprintf(fid, "title\n");
    fprintf(fid, "\n");
    fprintf(fid, "0 %d\n",mult);
  endif

  ## print the body of the z-matrix
  nvars = 0;
  svars = {};
  dvars = [];
  xzmat = zeros(rows(izmat),3);
  for i = 1:rows(izmat)
    ii0 = izmat(i,1);
    fprintf(fid,"%s ",mol.atname{ii0});
    if (izmat(i,2) > 0)
      i1 = izmat(i,2);
      ii1 = izmat(i1,1);
      a = op_dist(mol.atxyz(:,ii1),mol.atxyz(:,ii0));
      if (!isempty(izvar) && izvar(i,2))
        nvars++;
        fprintf(fid,"%d b%d ",i1,nvars);
        svars{nvars} = sprintf("b%d",nvars);
        dvars(nvars) = a;
      else
        fprintf(fid,"%d %.10f ",i1,a);
      endif
      xzmat(i,1) = a;
    endif
    if (izmat(i,3) > 0)
      i2 = izmat(i,3);
      ii2 = izmat(i2,1);
      a = op_angle(mol.atxyz(:,ii2),mol.atxyz(:,ii1),mol.atxyz(:,ii0));
      if (!isempty(izvar) && izvar(i,3))
        nvars++;
        fprintf(fid,"%d a%d ",i2,nvars);
        svars{nvars} = sprintf("a%d",nvars);
        dvars(nvars) = a;
      else
        fprintf(fid,"%d %.10f ",i2,a);
      endif
      xzmat(i,2) = a;
    endif
    if (izmat(i,4) > 0)
      i3 = izmat(i,4);
      ii3 = izmat(i3,1);
      a = op_dihedral(mol.atxyz(:,ii3),mol.atxyz(:,ii2),mol.atxyz(:,ii1),mol.atxyz(:,ii0));
      if (!isempty(izvar) && izvar(i,4))
        nvars++;
        fprintf(fid,"%d d%d ",i3,nvars);
        svars{nvars} = sprintf("d%d",nvars);
        dvars(nvars) = a;
      else
        fprintf(fid,"%d %.10f ",i3,a);
      endif
      xzmat(i,3) = a;
    endif
    fprintf(fid,"\n");
  endfor
  fprintf(fid,"\n");

  ## print the variables
  for i = 1:nvars
    fprintf(fid,"%s %.10f\n",svars{i},dvars(i));
  endfor
  fprintf(fid,"\n");

  ## close the file
  fclose(fid);

endfunction

