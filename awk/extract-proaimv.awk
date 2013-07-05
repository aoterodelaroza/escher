#! /usr/bin/awk -f
#
# extract-proaimv - Get atomic properties from the proaimv output file.
# This analyzer is ready to read many files and prepare a table of several
# properties per atom. Currently:
#
#     charges, dipole module and quadrupole eigenvalues
#
#
#-----------------------------------------------------------------------
# CopyRight (c): Victor Lua~na, Jul. 31, 1999
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

function bubble_sort(xvar, n,    i, j, xdummy) {
   for (i=1; i<=n-1; i++) {
      for (j=i+1; j<=n; j++) {
         if (xvar[i] > xvar[j]) {
            xdummy = xvar[i]
            xvar[i] = xvar[j]
            xvar[j] = xdummy
            }
         }
      }
   }

BEGIN {
   printf "# Analysis of PROAIMV output files\n#\n"
   printf "%-31s %-3s %3s", "# File", "Atm", "Num"
   printf " %12s %12s", "charge", "dipmod"
   printf " %12s", "Det.quad"
   printf " %12s %12s %12s", "quad1", "quad2", "quad3"
   printf " %12s\n", "eatom"

   printf "%-31s %-3s %3s", "# -----------------------------", "---", "---"
   printf " %12s %12s", "------------", "------------"
   printf " %12s", "------------"
   printf " %12s %12s %12s", "------------", "------------", "------------"
   printf " %12s\n", "------------"
   }

/^ PROAIMV/ { fich = FILENAME }

/INTEGRATION IS OVER ATOM/ { atom = $5; n_at = $6 }

/NET CHARGE/ { charge = $5 }

/E\(ATOM\)/ { eatom = $4 }

/ EL DX / { dip1 = $3 }
/ EL DY / { dip2 = $3 }
/ EL DZ / { dip3 = $3; dipmod = sqrt( dip1*dip1 + dip2*dip2 + dip3*dip3 ) }

/ EIGENVALUES OF QUADRUPOLE MOMENT TENSOR/ {
   getline
   quad[1] = $1
   quad[2] = $2
   quad[3] = $3
   bubble_sort(quad, 3)
   }

/NORMAL TERMINATION OF PROAIM/ {
   printf "%-31s %-3s %3d", fich, atom, n_at
   printf " %12.6f %12.6f", charge, dipmod
   printf " %12.6f", quad[1]*quad[2]*quad[3]
   printf " %12.6f %12.6f %12.6f", quad[1], quad[2], quad[3]
   printf " %12.6f\n", eatom
   nextfile
   }
