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

function err = cr_write_cif(cr, ciffile="default.cif", LOG=1)
% function err = cr_write_cif(cr, ciffile, LOG=1)
%
% cr_write_cif - writes the description of a crystal to a cif file.
%
% Required input variables:
% {cr}: crystal data. See cr_popinfo.m for a description of the data
% structure.
% {cif}: struct with data to be included in the final output. Required
% fields include:
%   cif.file ... cif output file to be created.
%   cif.code ... codename for this file.
% Other fields are optional and include:
%
% Optional input variables (all have default values):
% {ciffile}:
% {LOG}: print the final result if LOG>0.
%
% Ouput variables:
% {err}: error core (0: no error).
%
% 2012 Nov: VLC, first version
%

   bohrtoans = 0.52917720859;


   [fid,msg] = fopen(ciffile, 'w+');
   if (fid < 0 || ferror(fid))
      disp(msg)
      error("cr_write_cif: Could not open -- %s", ciffile);
   endif

   # cif headlines:
   code="default";
   fprintf(fid,"data_%s\n", code);
   fprintf(fid,";\n");

   # crystal physical properties:
   fprintf(fid,"_cell_volume %.6f\n", cr.omega*bohrtoans**3);
   fprintf(fid,";\n");

   # no symmetry is used in this starting version:
   fprintf(fid,"_symmetry_space_group_name_H-M 'P 1'\n");
   fprintf(fid,"symmetry_Int_Tables_number 1\n");
   fprintf(fid,"loop_\n");
   fprintf(fid,"_symmetry_equiv_pos_site_id\n");
   fprintf(fid,"_symmetry_equiv_pos_as_xyz\n");
   fprintf(fid,"1 x,y,z\n");
   fprintf(fid,"_cell_length_a %.9f\n", cr.a(1)*bohrtoans);
   fprintf(fid,"_cell_length_b %.9f\n", cr.a(2)*bohrtoans);
   fprintf(fid,"_cell_length_c %.9f\n", cr.a(3)*bohrtoans);
   fprintf(fid,"_cell_angle_alpha %.4f\n", cr.b(1)*180/pi);
   fprintf(fid,"_cell_angle_beta %.4f\n",  cr.b(2)*180/pi);
   fprintf(fid,"_cell_angle_gamma %.4f\n", cr.b(3)*180/pi);
   fprintf(fid,"_cell_formula_units_Z 1\n");
   fprintf(fid,"loop_\n");
   fprintf(fid,"_atom_type_symbol\n");
   fprintf(fid,"_atom_type_radius_bond\n");
   for i = 1:cr.ntyp
      zi=cr.ztyp(i);
      symb=cell2mat(cr.attyp(i));
      [Zat,atom]=mol_dbatom(symb);
      ri=atom.rcov;
      fprintf(fid,"%s %.3f\n", symb, ri);
   endfor
   fprintf(fid,"loop_\n");
   fprintf(fid,"_atom_site_label\n");
   fprintf(fid,"_atom_site_type_symbol\n");
   fprintf(fid,"_atom_site_fract_x\n");
   fprintf(fid,"_atom_site_fract_y\n");
   fprintf(fid,"_atom_site_fract_z\n");
   for i = 1:cr.nat
      ti=cr.typ(i);
      zi=cr.ztyp(ti);
      symb=strrep(cr.attyp{ti}," ","_");
      si=sprintf("%s%0d", symb, i);
      fprintf(fid,"%s %s %.9f %.9f %.9f\n", si, symb, cr.x(i,:));
   endfor
   fprintf(fid,"#END\n");
   err = 0;
   fclose(fid);

endfunction
