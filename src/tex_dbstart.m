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

function tex_dbstart();
% function tex_dbstart();
%
% tex_dbstart - start the internal texture database.

  global texdb
  
  texdb = cell(); n = 0;
  
  #### povray textures ####
  ## a simple default
  n = 1;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "tdefault";
  texdb{n}.string = "finish {phong 0.2 diffuse 0.5 ambient 0.3 reflection 0.03 brilliance 1.0}";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,%.4f,%.4f>}";

  ## default for balls, sticks and opaque triangles
  n = 2; 
  texdb{n} = texdb{1}; 
  texdb{n}.name = "ball_default";
  n = 3; 
  texdb{n} = texdb{1}; 
  texdb{n}.name = "stick_default";
  n = 4; 
  texdb{n} = texdb{1}; 
  texdb{n}.name = "opaque_triangle_default";

  ## povray textures, other
  n = 5;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "old_default";
  texdb{n}.string = "finish{specular 0.2 roughness 0.1 reflection 0.02}";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,%.4f,%.4f>}";

  n = 6;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "fglass5";
  texdb{n}.string = "finish { specular 0.7 roughness 0.001 ambient 0 diffuse 0 reflection {0.2, 1.0 fresnel on } conserve_energy }";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,0.8,0.0>}";

  n = 7;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "fglass10";
  texdb{n}.string = "finish { specular 0.6 roughness 0.002 ambient 0 diffuse 0.1 reflection {0.05, 1.0} conserve_energy }";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,0.8,0.0>}";

  n = 8;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "metal1";
  texdb{n}.string = "finish { brilliance 2 diffuse 0.391 ambient rgb <0.22000,0.18100,0.12100> reflection rgb <0.55000,0.45250,0.30250> metallic 1 specular 0.20 roughness 1/20}";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,%.4f,%.4f>}";

  n = 9;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "metal2";
  texdb{n}.string = "finish { ambient 0.35 brilliance 2 diffuse 0.3 metallic specular 0.80 roughness 1/20 reflection 0.1}";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,%.4f,%.4f>}";

  n = 10;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "metal3";
  texdb{n}.string = "finish { ambient 0.15 brilliance 5 diffuse 0.6 metallic specular 0.80 roughness 1/100 reflection 0.65}";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,%.4f,%.4f>}";

  n = 11;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "starfield";
  texdb{n}.string = "finish { diffuse 0 ambient 1 }";
  texdb{n}.pigment = "pigment { granite color_map { [ 0.000  0.260 color rgb < 0, 0, 0> color rgb < 0, 0, 0> ] [ 0.260  0.300 color rgb <.5,.5,.4> color rgb <.8,.8,.4> ] [ 0.300  0.460 color rgb < 0, 0, 0> color rgb < 0, 0, 0> ] [ 0.460  0.500 color rgb <.4,.4,.5> color rgb <.4,.4,.8> ] [ 0.500  0.660 color rgb < 0, 0, 0> color rgb < 0, 0, 0> ] [ 0.660  0.700 color rgb <.5,.4,.4> color rgb <.8,.4,.4> ] [ 0.700  0.860 color rgb < 0, 0, 0> color rgb < 0, 0, 0> ] [ 0.860  0.900 color rgb <.5,.5,.5> color rgb < 1, 1, 1> ] [ 0.900  1.000 color rgb < 0, 0, 0> color rgb < 0, 0, 0> ]} turbulence 1 sine_wave scale .5 }";

  n = 12;
  texdb{n} = texture();
  texdb{n}.typ = "pov";
  texdb{n}.name = "chrome";
  texdb{n}.string = "finish {ambient 0.3 diffuse 0.7 reflection 0.15 brilliance 8 specular 0.8 roughness 0.1}";
  texdb{n}.pigment = "pigment {color rgbft <%.4f,%.4f,%.4f,%.4f,%.4f>}";

  #### obj textures ####
  ## a simple default
  nstart = n;
  n = nstart + 1;
  texdb{n} = texture();
  texdb{n}.typ = "obj";
  texdb{n}.name = "tdefault";
  texdb{n}.Ns = 96.078;
  texdb{n}.Ka = [0 0 0];
  texdb{n}.Ks = [128 128 128];
  texdb{n}.Ni = 1;
  texdb{n}.illum = 2;

  ## default for balls, sticks and opaque triangles
  n = nstart + 2; 
  texdb{n} = texdb{nstart+1}; 
  texdb{n}.name = "ball_default";
  n = nstart + 3; 
  texdb{n} = texdb{nstart+2}; 
  texdb{n}.name = "stick_default";
  n = nstart + 4; 
  texdb{n} = texdb{nstart+3}; 
  texdb{n}.name = "opaque_triangle_default";

endfunction
