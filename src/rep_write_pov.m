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

function rep_write_pov(rep,file="",LOG=0)
% function rep_write_pov(rep,file,LOG=0)
%
% rep_write_pov - write a representation to a povray (pov) file.
%
% Required input variables:
% rep: representation.
% file: pov file. If no file is given, use standard output.
%
% Optional input variables (all have default values):
% {LOG}: verbose level (0=silent,1=verbose).
%

  bohrtoans = 0.52917720859;

  ## open file
  if (!isstruct(rep) || isempty(rep))
    error("Invalid rep input")
  endif
  if (!isempty(file))
    fid = fopen(file,"w");
  else
    fid = stdout();
  endif

  ## pov header
##  fprintf(fid,"#include \"colors.inc\"\n");
##  fprintf(fid,"#include \"glass.inc\"\n");
##  fprintf(fid,"#include \"textures.inc\"\n");
##  fprintf(fid,"#include \"woods.inc\"\n");
##  fprintf(fid,"#include \"stones.inc\"\n");
  fprintf(fid,"// Model created MolWare, adapted from tessel's writepov.f \n");
  if (isfield(rep,"name") && !isempty(rep.name))
    fprintf(fid,"// Name: %s \n",rep.name);
  endif
  fprintf(fid,"// \n\n");

  ## textures
  pigment = {};
  for i = 1:length(rep.texlib)
    tex = texture("pov",rep.texlib{i});
    pigment{i} = tex.pigment;
    fprintf(fid,"#declare %s = texture {%s}\n",rep.texlib{i},tex.string);
  endfor

  ## write balls to the pov file
  fprintf(fid,"#declare Mol1 = union {\n")
  for i = 1:rep.nball
    str = pigment{rep.ball{i}.tex};
    s = sprintf("%s %s %s","  sphere{<%.9f,%.9f,%.9f>, %.9f texture {%s",str,"}}\n");
    n = sum(str == "%");
    rgb = fillrgb(rep.ball{i}.rgb) / 255;
    rgb = rgb(1:n);
    if (!isempty(rgb))
      fprintf(fid,s,rep.ball{i}.x,rep.ball{i}.r,rep.texlib{rep.ball{i}.tex},rgb);
    else
      fprintf(fid,s,rep.ball{i}.x,rep.ball{i}.r,rep.texlib{rep.ball{i}.tex});
    endif
  endfor

  ## write sticks to the pov file
  for i = 1:rep.nstick
    ## skip degenerate cylinders
    if (norm(rep.stick{i}.x1 - rep.stick{i}.x0) < eps)
      continue
    endif
    str = pigment{rep.stick{i}.tex};
    s = sprintf("%s %s %s","   cylinder{<%.9f,%.9f,%.9f>,<%.9f,%.9f,%.9f>, %.9f texture {%s",str,"}}\n");
    n = sum(str == "%");
    rgb = fillrgb(rep.stick{i}.rgb) / 255;
    rgb = rgb(1:n);
    if (!isempty(rgb))
      fprintf(fid,s,rep.stick{i}.x0,rep.stick{i}.x1,rep.stick{i}.r,rep.texlib{rep.stick{i}.tex},rgb);
    else
      fprintf(fid,s,rep.stick{i}.x0,rep.stick{i}.x1,rep.stick{i}.r,rep.texlib{rep.stick{i}.tex});
    endif
  endfor

  for i = 1:rep.ntriangle
    str = pigment{rep.triangle{i}.tex};
    s = sprintf("%s %s %s","   triangle{<%.9f,%.9f,%.9f> <%.9f,%.9f,%.9f> <%.9f,%.9f,%.9f> texture {%s",str,"}}\n");
    n = sum(str == "%");
    rgb = fillrgb((rep.vertex{rep.triangle{i}.idx(1)}.rgb+...
                   rep.vertex{rep.triangle{i}.idx(2)}.rgb+...
                   rep.vertex{rep.triangle{i}.idx(3)}.rgb)/3)/255;
    rgb = rgb(1:n);
    if (!isempty(rgb))
      fprintf(fid,s,
              rep.vertex{rep.triangle{i}.idx(1)}.x,...
              rep.vertex{rep.triangle{i}.idx(2)}.x,...
              rep.vertex{rep.triangle{i}.idx(3)}.x,...
              rep.texlib{rep.triangle{i}.tex},rgb);
    else
      fprintf(fid,s,
              rep.vertex{rep.triangle{i}.idx(1)}.x,...
              rep.vertex{rep.triangle{i}.idx(2)}.x,...
              rep.vertex{rep.triangle{i}.idx(3)}.x,...
              rep.texlib{rep.triangle{i}.tex});
    end
  endfor

  ## end the global union
  fprintf(fid,"}\n\n");

  ## rotate
  fprintf(fid,"object {Mol1 rotate <0,0,0>}\n");
  fprintf(fid,"\n");

  ## surfaces
  for i = 1:rep.nsurf
    fprintf(fid,"mesh2 {\n");
    fprintf(fid,"  vertex_vectors {\n");
    fprintf(fid,"  %d,\n",size(rep.surf{i}.v,1));
    fprintf(fid,"  <%.10f,%.10f,%.10f>,\n",rep.surf{i}.v');
    fprintf(fid,"  }\n");
    fprintf(fid,"  normal_vectors {\n");
    fprintf(fid,"  %d,\n",size(rep.surf{i}.n,1));
    fprintf(fid,"  <%.10f,%.10f,%.10f>,\n",rep.surf{i}.n');
    fprintf(fid,"  }\n");
    fprintf(fid,"  texture_list {\n");
    fprintf(fid,"  1,\n");
    s = sprintf("  texture{%s %s}",rep.texlib{rep.surf{i}.ftex},pigment{rep.surf{i}.ftex});
    fprintf(fid,s,fillrgb(rep.surf{i}.frgb)/255);
    fprintf(fid,"\n  }\n");
    fprintf(fid,"  face_indices {\n");
    fprintf(fid,"  %d\n",size(rep.surf{i}.f,1));
    fprintf(fid,"  <%d,%d,%d>,0,\n",(rep.surf{i}.f-1)');
    fprintf(fid,"  }\n");
    fprintf(fid,"  inside_vector <0, 0, 1>\n");
    fprintf(fid,"}\n");
  endfor

  ## camera
  if (!isfield(rep.cam,"cop") || isempty(rep.cam.cop))
    error("empty camera: missing rep_setdefaultscene?");
  endif
  fprintf(fid,"camera {\n");
  if (isfield(rep.cam,"persp") && !isempty(rep.cam.persp))
    if (rep.cam.persp == 1)
      fprintf(fid,"  perspective\n");
    else
      fprintf(fid,"  orthographic\n");
    endif
  endif
  if (isfield(rep.cam,"cop") && !isempty(rep.cam.cop))
    fprintf(fid,"  location  <%.5f,%.5f,%.5f>\n",rep.cam.cop);
  endif
  if (isfield(rep.cam,"sky") && !isempty(rep.cam.sky))
    fprintf(fid,"  sky       <%.5f,%.5f,%.5f>\n",rep.cam.sky);
  endif
  if (isfield(rep.cam,"vuv") && !isempty(rep.cam.vuv))
    fprintf(fid,"  up        <%.5f,%.5f,%.5f>\n",rep.cam.vuv);
  endif
  if (isfield(rep.cam,"rht") && !isempty(rep.cam.rht))
    fprintf(fid,"  right     <%.5f,%.5f,%.5f>\n",rep.cam.rht);
  endif
  if (isfield(rep.cam,"drt") && !isempty(rep.cam.drt))
    fprintf(fid,"  direction <%.5f,%.5f,%.5f>\n",rep.cam.drt);
  endif
  if (isfield(rep.cam,"angle") && !isempty(rep.cam.angle))
    fprintf(fid,"  angle   %.5f\n",rep.cam.angle);
  endif
  if (isfield(rep.cam,"vrp") && !isempty(rep.cam.vrp))
    fprintf(fid,"  look_at   <%.5f,%.5f,%.5f>\n",rep.cam.vrp);
  endif
  if (isfield(rep.cam,"matrix") && !isempty(rep.cam.matrix))
    ## I use the COP as the initial position of the camera, then translate
    ## away using the modelview matrix. But for this to work I need to translate the
    ## COP to the origin before applying the translation.
    fprintf(fid,"  transform{\n");
    fprintf(fid,"    translate <%.10f, %.10f, %.10f>\n",rep.cam.cop);
    fprintf(fid,"  inverse }\n");
    fprintf(fid,"  transform{\n");
    fprintf(fid,"    matrix <\n");
    fprintf(fid,"    %.10f, %.10f, %.10f,\n",rep.cam.matrix(1,1:3));
    fprintf(fid,"    %.10f, %.10f, %.10f,\n",rep.cam.matrix(2,1:3));
    fprintf(fid,"    %.10f, %.10f, %.10f,\n",rep.cam.matrix(3,1:3));
    fprintf(fid,"    %.10f, %.10f, %.10f\n",rep.cam.matrix(4,1:3));
    fprintf(fid,"    >\n");
    fprintf(fid,"  inverse }\n");
    fprintf(fid,"  transform{\n");
    fprintf(fid,"    translate <%.10f, %.10f, %.10f>\n",rep.cam.cop);
    fprintf(fid,"  }\n");
  endif
  if (isfield(rep.cam,"rot") && !isempty(rep.cam.rot))
    fprintf(fid,"  rotate   <%.5f,%.5f,%.5f>\n",rep.cam.rot);
  endif
  if (isfield(rep.cam,"trans") && !isempty(rep.cam.trans))
    fprintf(fid,"  translate   <%.5f,%.5f,%.5f>\n",rep.cam.trans);
  endif
  fprintf(fid,"}\n");
  fprintf(fid,"\n");

  ## lights
  for i = 1:rep.nlight
    fprintf(fid,"light_source {\n  <%.5f,%.5f,%.5f>\n",rep.light{i}.x);
    fprintf(fid,"  color rgb<%.5f,%.5f,%.5f>\n",rep.light{i}.color/255*rep.light{i}.intensity);
    if (rep.light{i}.shadowless != 0)
      fprintf(fid,"  shadowless\n");
    endif
    if (isfield(rep.light{i},"matrix") && ismatrix(rep.light{i}.matrix) && !isempty(rep.light{i}.matrix))
      fprintf(fid,"  transform{\n");
      fprintf(fid,"    translate <%.10f, %.10f, %.10f>\n",rep.cam.cop);
      fprintf(fid,"  inverse }\n");
      fprintf(fid,"  transform{\n");
      fprintf(fid,"    matrix <\n");
      fprintf(fid,"    %.10f, %.10f, %.10f,\n",rep.light{i}.matrix(1,1:3));
      fprintf(fid,"    %.10f, %.10f, %.10f,\n",rep.light{i}.matrix(2,1:3));
      fprintf(fid,"    %.10f, %.10f, %.10f,\n",rep.light{i}.matrix(3,1:3));
      fprintf(fid,"    %.10f, %.10f, %.10f\n",rep.light{i}.matrix(4,1:3));
      fprintf(fid,"    >\n");
      fprintf(fid,"  inverse }\n");
      fprintf(fid,"  transform{\n");
      fprintf(fid,"    translate <%.10f, %.10f, %.10f>\n",rep.cam.cop);
      fprintf(fid,"  }\n");
    endif
    fprintf(fid,"}\n");
  endfor
  fprintf(fid,"\n");

  ## background color
  fprintf(fid," background {color rgb <%.5f,%.5f,%.5f>}\n",rep.bgcolor(1:3)/255);
  fprintf(fid,"\n");

  ## wrap up
  if (!isempty(file))
    root = strrep(file,".pov","");
    fprintf(fid,"//runme: povray -D -UV +I%s +O%s.png +W1000 +H1000 +A\n",file,root);
    fclose(fid);
  endif
  if (LOG > 0)
    printf("rep_write_pov: Writing %s\n", file);
  endif

endfunction
