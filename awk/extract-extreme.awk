#! /usr/bin/awk -f
#
# extract-extreme.awk - Extract bond critical point properties from
# the 'slightly verbose' output of extreme.
#
#-----------------------------------------------------------------------
# CopyRight (c): Victor Lua~na, Jul. 16, 1999
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
#

BEGIN {
   eldens = 0
   type_deg = -50
   type_nuc = -40
   type_bnd = -30
   type_rng = -20
   type_min = -10
   type_nnm = -60
   critpoint[type_deg] = "Deg?"
   critpoint[type_nuc] = "Nuc."
   critpoint[type_bnd] = "Bond"
   critpoint[type_rng] = "Ring"
   critpoint[type_min] = "Cage"
   critpoint[type_nnm] = "NNM "
   N_deg = 0
   N_nuc = 0
   N_bnd = 0
   N_rng = 0
   N_min = 0
   N_nnm = 0
   nother = 5
   typother[1] = type_nuc
   typother[2] = type_rng
   typother[3] = type_min
   typother[4] = type_nnm
   typother[5] = type_deg
   }

/^ EXTREME/ {
   getline
   getline
   runtitle = $0
   }

/^ NUCLEAR COORDINATES/,/^ *$/ {
   if (NF==5) {
      nat++
      if (nat!=$2) {print "Atom number error!!!"}  # exit[1]
      atname[nat] = $1
      atnumber[nat] = $2
      atx[nat] = $3
      aty[nat] = $4
      atz[nat] = $5
      }
   }

/^ ANALYZING RHO, THE ELECTRON DENSITY/ { eldens = 1 }

/^ NEW CRITICAL POINT FOUND/ {
   ncrit++
   nneig = 0
   }

/^ COORDINATES OF CRITICAL POINT AND DISTANCE FROM MOLECULAR ORIGIN/ {
   getline; xcrit[ncrit] = $3
   getline; ycrit[ncrit] = $3
   getline; zcrit[ncrit] = $3
   }

/^ VECTORS AND DISTANCES FROM NUCLEI TO CRITICAL POINT/,/^ ANGLES/ {
   if (NF==6) {
      nneig++
      if (nneig!=$2) {print "Atom number error 2!!!"}  #  exit[1]
      dist[nneig] = $6
      iord[nneig] = nneig
      }
   }

/^ EIGENVALUES OF THE HESSIAN/ {
   getline
   hess1[ncrit] = $1
   hess2[ncrit] = $2
   hess3[ncrit] = $3
   #
   # Classify the critical point:
   #
   rank = 0
   signature = 0
   for (i=1; i<=3; i++) {
      if ($i > 0) { rank++; signature++ }
      else if ($i < 0) { rank++; signature-- }
      }
   if (rank < 3) { type[ncrit] = type_deg; N_deg++ }
   else {
      if (signature == -3) { type[ncrit] = type_nuc; N_nuc++ }
      else if (signature == -1) { type[ncrit] = type_bnd; N_bnd++ }
      else if (signature == +1) { type[ncrit] = type_rng; N_rng++ }
      else if (signature == +3) { type[ncrit] = type_min; N_min++ }
      }
   crittype[ncrit] = sprintf("(%d,%d)", rank, signature)
   #
   # If this is a bond cp sort the neighbors and guess the atoms bonded:
   #
   if (type[ncrit] == type_bnd) {
      for (i=1; i<=nneig; i++) {
         for (j=1; j<i; j++) {
            if (dist[iord[j]] > dist[iord[i]]) {
               k = iord[i]
               iord[i] = iord[j]
               iord[j] = k
               }
            }
         }
      neig1[ncrit] = iord[1]
      neig2[ncrit] = iord[2]
      dist1[ncrit] = dist[iord[1]]
      dist2[ncrit] = dist[iord[2]]
      }
   #
   # If this is a maximum guess whether it is a nucleus or a nnm
   # The criterion is just a distance: the nuclear maximum is allowed
   # to be up to 0.05 bohr apart from the nucleus. This separation
   # occurs typically in hydrogens.
   #
   if (type[ncrit] == type_nuc) {
      for (i=1; i<=nneig; i++) {
         for (j=1; j<i; j++) {
            if (dist[iord[j]] > dist[iord[i]]) {
               k = iord[i]
               iord[i] = iord[j]
               iord[j] = k
               }
            }
         }
      neig1[ncrit] = iord[1]
      dist1[ncrit] = dist[iord[1]]
      if (dist1[ncrit] > 5e-2) { type[ncrit] = type_nnm; N_nnm++ }
      if (atnumber[neig1[ncrit]] < 3 && dist1[ncrit] < 1e-1) {
         type[ncrit] = type_nuc
         }
      }
   }

/^ Rho\(r\)/             { rho[ncrit] = $2 }
/^ DEL\*\*2\(Rho\(r\)\)/ { del2[ncrit] = $2 }
/^ G\(r\)/               { kinG[ncrit] = $2 }

END {
   #
   # Safety check: Does this file correspond to the topology of
   # the electron density?
   #
   if (eldens != 1) {
      print "WARNING!"
      print "WARNING: This appears not to correspond to the Electron Density"
      print "WARNING!"
      }

   printf "Data extracted from the file: %s\n\n", FILENAME
   printf "%s\n\n", runtitle
   printf "Atoms and positions\n"
   for (i=1; i<=nat; i++) {
      printf "%3d ", i
      printf "%-6s ", atname[i]
      printf "%12.6f ", atx[i]
      printf "%12.6f ", aty[i]
      printf "%12.6f ", atz[i]
      printf "\n"
      }
   printf "\n"

   #
   # Sort the critical points from high to low density values
   #
   for (i=1; i<=ncrit; i++) { isort[i] = i }
   for (i=1; i<ncrit; i++) {
      for (j=i+1; j<=ncrit; j++) {
         if (rho[isort[i]] < rho[isort[j]]) {
            idum = isort[j]
            isort[j] = isort[i]
            isort[i] = idum
            }
         }
      }

   printf "Bond critical points\n"
   printf "Num "
   printf "Atom     "
   printf "Atom     "
   printf "    Radius 1 "
   printf "    Radius 2 "
   printf "       Critical Point position (x,y,z)"
   printf "   E. density"
   printf "    Laplacian"
   printf "    G: Kin.E."
   printf "          Hessian matrix eigenvalues  "
   printf "\n"
   for (i=1; i<=ncrit; i++) {
      i1 = isort[i]
      if (type[i1] == type_bnd) {
         printf "%3d ", i1
         printf "%-4s(%02d) ", atname[neig1[i1]], neig1[i1]
         printf "%-4s(%02d) ", atname[neig2[i1]], neig2[i1]
         printf "%12.6f %12.6f ", dist1[i1], dist2[i1]
         printf "%12.6f %12.6f %12.6f ", xcrit[i1], ycrit[i1], zcrit[i1]
         printf "%12.6f %12.6f %12.6f ", rho[i1], del2[i1], kinG[i1]
         printf "%12.6f %12.6f %12.6f ", hess1[i1], hess2[i1], hess3[i1]
         printf "\n"
         }
      }
   printf "\n"

   printf "Other critical points\n"
   printf "Num "
   printf "Type   "
   printf "Atom     "
   printf "(r,s)  "
   printf "       Critical Point position (x,y,z)"
   printf "   E. density"
   printf "    Laplacian"
   printf "    G: Kin.E."
   printf "          Hessian matrix eigenvalues  "
   printf "\n"
   for (k=1; k<=nother; k++) {
      for (i=1; i<=ncrit; i++) {
         i1 = isort[i]
         if (type[i1] == typother[k]) {
            printf "%3d ", i1
            printf "%-6s ", critpoint[type[i1]]
            if (type[i1] == type_nuc) {
               printf "%-4s(%02d) ", atname[neig1[i1]], neig1[i1]
               }
            else { printf "         " }
            printf "%-6s ", crittype[i1]
            printf "%12.6f %12.6f %12.6f ", xcrit[i1], ycrit[i1], zcrit[i1]
            printf "%12.6f ", rho[i1]
            if (type[i1] != type_nuc) {
               printf "%12.6f %12.6f ", del2[i1], kinG[i1]
               printf "%12.6f %12.6f %12.6f ", hess1[i1], hess2[i1], hess3[i1]
               }
            printf "\n"
            }
         }
      }
   printf "\n"

   # Printing the neighborhood of NNM:
   if (N_nnm > 0) {
      printf "\n"
      for (i=1; i<=ncrit; i++) {
         i1 = isort[i]
         if (type[i1] == type_nnm) {
            printf "Neighbors for NNM %4d: ", i1
            for (j=1; j<=nat; j++) {
               idis[j] = j
               xx = atx[j] - xcrit[i1]
               yy = aty[j] - ycrit[i1]
               zz = atz[j] - zcrit[i1]
               dis[j] = sqrt( xx*xx + yy*yy + zz*zz )
               }
            # Sort the neighbors by distances:
            for (k=1; k<=nat; k++) {
               for (l=k+1; l<=nat; l++) {
                  if (dis[idis[k]] > dis[idis[l]]) {
                     idum = idis[k]
                     idis[k] = idis[l]
                     idis[l] = idum
                     }
                  }
               }
            # Print out the first four neighbors:
            for (k=1; k<=4; k++) {
               k1 = idis[k]
               printf "%10.4f%4s%4d", dis[k1], atname[k1], k1
               }
            printf "\n"
            }
         }
      printf "\n"
      }

   printf "Total number of c.p.: %4d\n", ncrit
   printf "Types of c.p. (nbrc): %4d %4d %4d %4d\n", N_nuc, N_bnd, N_rng, N_min
   printf "Non-nuclear maxima  : %4d\n", N_nnm
   printf "Degenerate c.p.     : %4d\n", N_deg
   }
