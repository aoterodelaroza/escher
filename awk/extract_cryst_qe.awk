#! /usr/bin/awk -f

#
# extract_cryst_qe.awk - Analysis of the quantum espresso output that
# corresponds to the optimization of a molecular crystal. The last geometry
# is converted to a tessel input.
#
#-----------------------------------------------------------------------
# CopyRight (c): Victor Lua~na, Jun. 2011
#                Departamento de Quimica Fisica y Analitica
#                Universidad de Oviedo, 33007-Oviedo, Spain
#                victor@carbono.quimica.uniovi.es
#
# This script is distributed as an Open Source code protected by the
# GNU General Public License version 2 (GPL2) as published by the
# Free Software Foundation (http://www.gnu.org).
# In short: you can make and distribute copies of this script; you can
# also make changes in the code and distribute the modified files; but
# you should always document the introduced changes and you must distribute
# or make available to others the source codes.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#-----------------------------------------------------------------------
#

BEGIN { pi = atan2(1,1)*4 }

function myabs(x) {
   if (x < 0.0) { return -x }
   else { return x }
   endif
}
function myasin(x) { return atan2(x, sqrt(1-x*x)) }
function myacos(x) { return atan2(sqrt(1-x*x), x) }
function myatan(x) { return atan2(x,1) }

FNR==1 { doit = 0 }

/Begin final coordinates/ { doit = 1 }

/CELL_PARAMETERS/ {
   if (doit==1) {
      getline; r[1,1]=$1; r[1,2]=$2; r[1,3]=$3;
      getline; r[2,1]=$1; r[2,2]=$2; r[2,3]=$3;
      getline; r[3,1]=$1; r[3,2]=$2; r[3,3]=$3;
      nat = 0
      for (i=1; i<=3; i++) {
         for (j=1; j<=3; j++) {
            g[i,j] = 0
            for (k=1; k<=3; k++) {
               g[i,j] += r[i,k] * r[j,k]
            }
         }
      }
      aa = sqrt(g[1,1])
      bb = sqrt(g[2,2])
      cc = sqrt(g[3,3])
      alp = myacos(g[1,2] / (aa*bb)) * 180/pi
      bet = myacos(g[1,3] / (aa*cc)) * 180/pi
      gam = myacos(g[2,3] / (bb*cc)) * 180/pi
      print pi, aa, bb, cc, alp, bet, gam
   }
}

/ATOMIC_POSITIONS/,/End final coordinates/ {
   if (doit==1 && NF==4 ) {
      at_name[++nat] = $1
      at_x[nat] = $2
      at_y[nat] = $3
      at_z[nat] = $4
      print nat
   }
}

/lattice parameter (alat)/ {
   alat = $(NF-1)
   if (myabs(alat-1) > 1e-4) {
      printf "ERROR: alat must be 1 for this to work correctly!"
      exit(1)
   }
}

###/JOB DONE/ {}

END {
fi = "cr.tess";
printf "set camera COP -8 0 22 SKY 0 0 1 VUV 0 0 1 RHT 1 0 0 DRT 0 10 0\n" > fi
printf "set use_planes .false.\n" >> fi
printf "set background background {color rgb <1.0,1.0,1.0>}\n" >> fi
printf "set LENGTHUNIT angstrom\n" >> fi
printf "molecule\n" >> fi
printf "crystal\n" >> fi
printf "spg p 1\n" >> fi
printf "cell %.6f %.6f %.6f %.6f %.6f %.6f\n", aa, bb, cc, alp, bet, gam >> fi
for (i=1; i<=nat; i++) {
   printf "neq %.6f %.6f %.6f %s # %d\n", \
          at_x[i], at_y[i], at_z[i], at_name[i], i >> fi
}
printf "crystalbox  -1.01 -1.01 -1.30 2.01 2.01 2.30\n" >> fi
printf "clippingbox -0.01 -0.01 -0.01 1.01 1.01 1.01\n" >> fi
printf "endcrystal\n" >> fi
printf "set ball_texture Dull\n" >> fi
printf "set stick_texture Dull\n" >> fi
printf "set glass_texture pigment {color Clear} &\n" >> fi
printf "    finish {ambient 0.2 diffuse 0.2 specular 0.0 &\n" >> fi
printf "    reflection 0.0 roughness 0.02 refraction 1 ior 1.005}\n" >> fi
printf "set pov_texture 3 pigment {color Clear} finish {ambient 0.2 &\n" >> fi
printf "    diffuse 0.2 specular 0.4 reflection 0.1 roughness 0.02 & \n" >> fi
printf "    refraction 1 ior 1.005} normal {bumps 0.5 scale 0.05}\n" >> fi
printf "unitcell radius 0.1 rgb 1.0 0.5 0.5\n" >> fi
printf "set atomicdata F cpk 15\n" >> fi
printf "molmotif allmaincell jmol\n" >> fi
printf "povray cr.pov\n" >> fi
printf "off cr.off\n" >> fi
printf "vrml cr.wrl\n" >> fi
printf "endmolecule\n" >> fi
printf "run povray -D +Icr.pov +Ocr.png +W900 +H900 +A\n" >> fi
printf "end\n" >> fi
}
