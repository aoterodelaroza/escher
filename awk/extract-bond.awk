#! /usr/bin/awk -f
#
# extract-bond.awk - Use the S_ij(\omega) overlap matrix from the PROAIMV
# output to determine the number of electron pairs localized into a center
# and the number of pairs involved in a bond.
#
#-----------------------------------------------------------------------
# CopyRight (c): Victor Lua~na, Feb. 23, 2000
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
# RUN:  extract-bond.awk  INPUT > OUTPUT
#
# Content of the INPUT file:
#
#   TITLE  title_for_the_run.s
#
#   ATOM   name.s   Z.i   ncopies.i    proaimv_file.s
#     Several lines describing the atomic composition of the molecule:
#     name .......... a name for the atom (used only in labelling)
#     Z ............. atomic number
#     ncopies ....... how many atoms equivalent to this occur in the molecule
#     proaimv_file .. name of the file containing the proaimv output
#
#-----------------------------------------------------------------------
#
# Input file example: "h2npcl4.bond" --->
#
# TITLE H2NPCl4 Opt HF/6-31G**
# ATOM  H   1  1  h2npcl4.pm-h1.int
# ATOM  H   1  1  h2npcl4.pm-h2.int
# ATOM  N   7  1  h2npcl4.pm-n3.int
# ATOM  P  15  1  h2npcl4.pm-p4.int
# ATOM  Cl 17  1  h2npcl4.pm-cl5.int
# ATOM  Cl 17  1  h2npcl4.pm-cl6.int
# ATOM  Cl 17  1  h2npcl4.pm-cl7.int
# ATOM  Cl 17  1  h2npcl4.pm-cl8.int
#
#-----------------------------------------------------------------------

function NumberOfOrbitals(n, ppp) {
   ppp = int(sqrt(n))
   while ((ppp*(ppp+1))/2 < n) { ppp++ }
   return (ppp)
   }

function abs(x) { return( x>=0 ? x : -x ) }

BEGIN {
   IGNORECASE = 1
   natom = 0
   allatoms = 0
   }

#
# Read in the general title:
#
$1 == "TITLE"  { $1 = ""; title = $0 }

#
# Read a new atom in the molecule and recover the topological data from
# its proaimv output file:
#
$1 == "ATOM"   {
   natom++
   atname[natom] = $2
   atZ[natom]    = $3
   atcopy[natom] = $4
   atfile[natom] = $5
   for (i=1; i<=atcopy[natom]; i++) {allatoms++; allIA[allatoms] = natom}
   #
   # Start reading the proaimv output file:
   #
   while (getline < atfile[natom]) {
      if ($1 == "N") {
         atNel[natom] = $2
         atCharge[natom] = $5
         }
      #
      # Read in the Overlap matrix:
      # Assumption: overlap among orbitals from a restricted wf
      #
      if (/The Atomic Overlap Matrix/ ) {
         getline < atfile[natom]
         getline < atfile[natom]
         getline < atfile[natom]
         getline < atfile[natom]
         ndat = 0
         while (NF>0) {
            for (i=1; i<=NF; i++) { ndat++; data[ndat]=$i }
            getline < atfile[natom]
            }
         atNorb[natom] = NumberOfOrbitals(ndat,pqrst)
         k = 0
         for (i=1; i<=atNorb[natom]; i++) {
            for (j=1; j<=i; j++) {
               k++
               atS[natom,i,j] = data[k]
               atS[natom,j,i] = data[k]
               }
            }
         }
      }
      close(atfile[natom])
#mdc*if DEBUG output:
#   print "ATOM:", natom, atname[natom], atfile[natom]
#   print "-- Electrons:", atNel[natom]
#   print "-- Overlap:"
#   for (i=1; i<=atNorb[natom]; i++) {
#      for (j=1; j<=i; j++) {
#         printf "%8.3f", atS[natom,i,j]
#         }
#      printf "\n"
#      }
#mdc*endif DEBUG
   }

#
# Analyze the data and output the results:
#
END {
   for (i1=1; i1<=allatoms; i1++) {
      i = allIA[i1]
      for (j1=1; j1<=i1; j1++) {
         j = allIA[j1]
         F[i1,j1] = 0
         for (k=1; k<=atNorb[i]; k++) {
            for (l=1; l<=atNorb[j]; l++) {
               F[i1,j1] -= atS[i,k,l] * atS[j,k,l]
               }
            }
         F[i1,j1] = 2*F[i1,j1]
         F[j1,i1] = F[i1,j1]
         }
      }
   for (i1=1; i1<=allatoms; i1++) {
      atlambda[i1] = abs(F[i1,i1])
      for (j1=1; j1<i1; j1++) {
         delta[i1,j1] = abs(F[i1,j1]) + abs(F[j1,i1])
         delta[j1,i1] = abs(F[i1,j1]) + abs(F[j1,i1])
         }
      }

   print "BOND ANALYSIS: "
   print "Determine the number of electron pairs localized in an atomic"
   print "basin and the number of pairs in each bond."
   print ""
   print "Pair localized on the atoms:"
   print "  i Atom           Z  copy    lambda(A)            N(A)    Norb Proaimv-file"
   for (i1=1; i1<=allatoms; i1++) {
      i = allIA[i1]
      printf "%3d %-10s%6d%6d", i, atname[i], atZ[i], atcopy[i]
      printf "%15.6f%15.6f", atlambda[i1], atNel[i]
      printf "%6d  %s\n", atNorb[i], atfile[i]
      }

   print ""
   print "Pairs localized into the bonds:"
   printf "     "
   for (i1=1; i1<=allatoms; i1++) {
      i = allIA[i1]
      printf "%15s", atname[i]
      }
   printf "\n"
   for (i1=2; i1<=allatoms; i1++) {
      i = allIA[i1]
      printf " %-4s", atname[i]
      for (j1=1; j1<i1; j1++) {
         j = allIA[j1]
         printf "%15.6f", delta[i1,j1]
         }
      printf "\n"
      }

   print ""
   print "D2 matrix:"
   for (i1=1; i1<=allatoms; i1++) {
      i = allIA[i1]
      for (j1=1; j1<=i1; j1++) {
         j = allIA[j1]
         printf "%15.6f", (atNel[i]*atNel[j]+F[i1,j1])/2
         }
      printf "\n"
      }

   print ""
   print "F(i,j) matrix:"
   for (i1=1; i1<=allatoms; i1++) {
      i = allIA[i1]
      for (j1=1; j1<=i1; j1++) {
         j = allIA[j1]
         printf "%15.6f", F[i1,j1]
         }
      printf "\n"
      }
   
   }
