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

function PrintAngle(i1, i2, i3){
   x21 = atxc[i2] - atxc[i1]
   y21 = atyc[i2] - atyc[i1]
   z21 = atzc[i2] - atzc[i1]
   r21 = sqrt(x21*x21 + y21*y21 + z21*z21)
   x23 = atxc[i2] - atxc[i3]
   y23 = atyc[i2] - atyc[i3]
   z23 = atzc[i2] - atzc[i3]
   r23 = sqrt(x23*x23 + y23*y23 + z23*z23)
   pesc = (x21*x23 + y21*y23 + z21*z23)
   x123 = y21*z23 - y23*z21
   y123 = z21*x23 - z23*x21
   z123 = x21*y23 - x23*y21
   pvec = sqrt(x123^2 + y123^2 + z123^2)
   angle = atan2(pvec,pesc)
   angle = 45 * angle / atan2(1,1)
   lbl0 = sprintf("%s-", atname[i1])
   lbl0 = lbl0 sprintf("%s-", atname[i2])
   lbl0 = lbl0 sprintf("%s", atname[i3])
   lbl1 = sprintf("%s(%d)-", atname[i1], i1)
   lbl1 = lbl1 sprintf("%s(%d)-", atname[i2], i2)
   lbl1 = lbl1 sprintf("%s(%d)", atname[i3], i3)
   printf "ANGLE %12.6f %-12s %s\n", angle, lbl0, lbl1
   }

BEGIN {
   bohr2angstrom = 0.5291771
   bondmax = 1.15
   atomicnumber["H" ] =   1;    rcov[  1] =  53;
   atomicnumber["HE"] =   2;    rcov[  2] =  70;
   atomicnumber["LI"] =   3;    rcov[  3] =  68;
   atomicnumber["BE"] =   4;    rcov[  4] =  35;
   atomicnumber["B" ] =   5;    rcov[  5] =  83;
   atomicnumber["C" ] =   6;    rcov[  6] =  68;
   atomicnumber["N" ] =   7;    rcov[  7] =  68;
   atomicnumber["O" ] =   8;    rcov[  8] =  68;
   atomicnumber["F" ] =   9;    rcov[  9] =  64;
   atomicnumber["NE"] =  10;    rcov[ 10] =  70;
   atomicnumber["NA"] =  11;    rcov[ 11] =  97;
   atomicnumber["MG"] =  12;    rcov[ 12] = 110;
   atomicnumber["AL"] =  13;    rcov[ 13] = 135;
   atomicnumber["SI"] =  14;    rcov[ 14] = 120;
   atomicnumber["P" ] =  15;    rcov[ 15] = 105;
   atomicnumber["S" ] =  16;    rcov[ 16] = 102;
   atomicnumber["CL"] =  17;    rcov[ 17] =  99;
   atomicnumber["AR"] =  18;    rcov[ 18] =  70;
   atomicnumber["K" ] =  19;    rcov[ 19] = 133;
   atomicnumber["CA"] =  20;    rcov[ 20] =  99;
   atomicnumber["SC"] =  21;    rcov[ 21] = 144;
   atomicnumber["TI"] =  22;    rcov[ 22] = 147;
   atomicnumber["V" ] =  23;    rcov[ 23] = 133;
   atomicnumber["CR"] =  24;    rcov[ 24] = 135;
   atomicnumber["MN"] =  25;    rcov[ 25] = 135;
   atomicnumber["FE"] =  26;    rcov[ 26] = 134;
   atomicnumber["CO"] =  27;    rcov[ 27] = 133;
   atomicnumber["NI"] =  28;    rcov[ 28] = 150;
   atomicnumber["CU"] =  29;    rcov[ 29] = 152;
   atomicnumber["ZN"] =  30;    rcov[ 30] = 145;
   atomicnumber["GA"] =  31;    rcov[ 31] = 122;
   atomicnumber["GE"] =  32;    rcov[ 32] = 117;
   atomicnumber["AS"] =  33;    rcov[ 33] = 121;
   atomicnumber["SE"] =  34;    rcov[ 34] = 122;
   atomicnumber["BR"] =  35;    rcov[ 35] = 121;
   atomicnumber["KR"] =  36;    rcov[ 36] = 191;
   atomicnumber["RB"] =  37;    rcov[ 37] = 147;
   atomicnumber["SR"] =  38;    rcov[ 38] = 112;
   atomicnumber["Y" ] =  39;    rcov[ 39] = 178;
   atomicnumber["ZR"] =  40;    rcov[ 40] = 157;
   atomicnumber["NB"] =  41;    rcov[ 41] = 148;
   atomicnumber["MO"] =  42;    rcov[ 42] = 147;
   atomicnumber["TC"] =  43;    rcov[ 43] = 135;
   atomicnumber["RU"] =  44;    rcov[ 44] = 140;
   atomicnumber["RH"] =  45;    rcov[ 45] = 145;
   atomicnumber["PD"] =  46;    rcov[ 46] = 150;
   atomicnumber["AG"] =  47;    rcov[ 47] = 159;
   atomicnumber["CD"] =  48;    rcov[ 48] = 169;
   atomicnumber["IN"] =  49;    rcov[ 49] = 163;
   atomicnumber["SN"] =  50;    rcov[ 50] = 146;
   atomicnumber["SB"] =  51;    rcov[ 51] = 146;
   atomicnumber["TE"] =  52;    rcov[ 52] = 147;
   atomicnumber["I" ] =  53;    rcov[ 53] = 140;
   atomicnumber["XE"] =  54;    rcov[ 54] = 198;
   atomicnumber["CS"] =  55;    rcov[ 55] = 167;
   atomicnumber["BA"] =  56;    rcov[ 56] = 134;
   atomicnumber["LA"] =  57;    rcov[ 57] = 187;
   atomicnumber["CE"] =  58;    rcov[ 58] = 183;
   atomicnumber["PR"] =  59;    rcov[ 59] = 182;
   atomicnumber["ND"] =  60;    rcov[ 60] = 181;
   atomicnumber["PM"] =  61;    rcov[ 61] = 180;
   atomicnumber["SM"] =  62;    rcov[ 62] = 180;
   atomicnumber["EU"] =  63;    rcov[ 63] = 199;
   atomicnumber["GD"] =  64;    rcov[ 64] = 179;
   atomicnumber["TB"] =  65;    rcov[ 65] = 176;
   atomicnumber["DY"] =  66;    rcov[ 66] = 175;
   atomicnumber["HO"] =  67;    rcov[ 67] = 174;
   atomicnumber["ER"] =  68;    rcov[ 68] = 173;
   atomicnumber["TM"] =  69;    rcov[ 69] = 172;
   atomicnumber["YB"] =  70;    rcov[ 70] = 194;
   atomicnumber["LU"] =  71;    rcov[ 71] = 172;
   atomicnumber["HF"] =  72;    rcov[ 72] = 157;
   atomicnumber["TA"] =  73;    rcov[ 73] = 143;
   atomicnumber["W" ] =  74;    rcov[ 74] = 137;
   atomicnumber["RE"] =  75;    rcov[ 75] = 135;
   atomicnumber["OS"] =  76;    rcov[ 76] = 137;
   atomicnumber["IR"] =  77;    rcov[ 77] = 132;
   atomicnumber["PT"] =  78;    rcov[ 78] = 150;
   atomicnumber["AU"] =  79;    rcov[ 79] = 150;
   atomicnumber["HG"] =  80;    rcov[ 80] = 170;
   atomicnumber["TL"] =  81;    rcov[ 81] = 155;
   atomicnumber["PB"] =  82;    rcov[ 82] = 154;
   atomicnumber["BI"] =  83;    rcov[ 83] = 154;
   atomicnumber["PO"] =  84;    rcov[ 84] = 168;
   atomicnumber["AT"] =  85;    rcov[ 85] = 170;
   atomicnumber["RN"] =  86;    rcov[ 86] = 240;
   atomicnumber["FR"] =  87;    rcov[ 87] = 200;
   atomicnumber["RA"] =  88;    rcov[ 88] = 190;
   atomicnumber["AC"] =  89;    rcov[ 89] = 188;
   atomicnumber["TH"] =  90;    rcov[ 90] = 179;
   atomicnumber["PA"] =  91;    rcov[ 91] = 161;
   atomicnumber["U" ] =  92;    rcov[ 92] = 158;
   atomicnumber["NP"] =  93;    rcov[ 93] = 155;
   atomicnumber["PU"] =  94;    rcov[ 94] = 153;
   atomicnumber["AM"] =  95;    rcov[ 95] = 151;
   atomicnumber["CM"] =  96;    rcov[ 96] = 151;
   atomicnumber["BK"] =  97;    rcov[ 97] = 151;
   atomicnumber["CF"] =  98;    rcov[ 98] = 151;
   atomicnumber["ES"] =  99;    rcov[ 99] = 151;
   atomicnumber["FM"] = 100;    rcov[100] = 151;
   atomicnumber["MD"] = 101;    rcov[101] = 151;
   atomicnumber["NO"] = 102;    rcov[102] = 151;
   atomicnumber["LW"] = 103;    rcov[103] = 151;
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
      atZ[nat] = atomicnumber[toupper($1)]
      atxc[nat] = $3
      atyc[nat] = $4
      atzc[nat] = $5
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
#   if (type[ncrit] == type_bnd) {
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
#      }
   #
   # If this is a maximum guess whether it is a nucleus or a nnm
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
      if (dist1[ncrit] > 1e-2) { type[ncrit] = type_nnm }
      }
   }

/^ Rho\(r\)/             { rho[ncrit] = $2 }
/^ DEL\*\*2\(Rho\(r\)\)/ { del2[ncrit] = $2 }
/^ G\(r\)/               { kinG[ncrit] = $2 }

END {
   printf "FICHERO %s\n\n", FILENAME
   printf "TITULO %s\n\n", runtitle
   printf "#Atoms and positions\n"
   for (i=1; i<=nat; i++) {
      printf "ATOM "
      printf "%-6s ", atname[i]
      printf "%12.6f ", atxc[i]
      printf "%12.6f ", atyc[i]
      printf "%12.6f ", atzc[i]
      printf "\n"
      }
   printf "\n"


   # Determine bonded atoms and print bond distances
   for (i=1; i<=nat; i++) {
      for (j=1; j<i; j++) {
         xx = atxc[i] - atxc[j]
         yy = atyc[i] - atyc[j]
         zz = atzc[i] - atzc[j]
         dij = sqrt(xx*xx + yy*yy + zz*zz) * bohr2angstrom
         dcov = (rcov[atZ[i]] + rcov[atZ[j]]) / 100
         bond[i,j] = (dij <= bondmax * dcov)
         bond[j,i] = bond[i,j]
         if (bond[i,j]) {
            i1 = i; i2 = j
            lbl0 = sprintf("%s-", atname[i1])
            lbl0 = lbl0 sprintf("%s", atname[i2])
            lbl1 = sprintf("%s(%d)-", atname[i1], i1)
            lbl1 = lbl1 sprintf("%s(%d)", atname[i2], i2)
            printf "BOND %12.6f %-9s %s\n", dij, lbl0, lbl1
            }
         }
      }
   printf "\n"

   # Print significative bond angles
   for (i=1; i<=nat; i++) {
      for (j=1; j<i; j++) {
         for (k=1; k<j; k++) {
            if (bond[i,j] && bond[j,k]) { PrintAngle(i,j,k) }
            if (bond[i,k] && bond[k,j]) { PrintAngle(i,k,j) }
            if (bond[j,i] && bond[i,k]) { PrintAngle(j,i,k) }
            }
         }
      }
   printf "\n"

   printf "#Bond critical points\n"
   printf "#       "
   printf "Atom     "
   printf "Atom     "
   printf " Radius 1 "
   printf " Radius 2 "
   printf "   E. density"
   printf "    Laplacian"
   printf "  Lambda.Perp"
   printf "  Epsilon"
   printf "    G: Kin.E."
   printf "\n"
   for (i=1; i<=ncrit; i++) {
      if (type[i] == type_bnd) {
         epsilon = 100 * (1 - hess2[i] / hess1[i])
         printf "CP BOND "
         printf "%-4s(%02d) ", atname[neig1[i]], neig1[i]
         printf "%-4s(%02d) ", atname[neig2[i]], neig2[i]
         printf "%9.3f %9.3f ", dist1[i], dist2[i]
         printf "%12.6f %12.6f ", rho[i], del2[i]
         printf "%12.6f %8.2f ",  hess3[i], epsilon
         printf "%12.6f ", kinG[i]
         printf "\n"
         #MDC#if debug
            printf "coord: %12.6f %12.6f %12.6f\n", xcrit[i], ycrit[i], zcrit[i]
         #MDC#endif debug
         }
      }
   printf "\n"

   printf "#Other critical points\n"
   printf "#         "
   printf "Type   "
   printf "Atom     "
   printf "(r,s)  "
   printf "   E. density"
   printf "    Laplacian"
   printf "  Lambda.Perp"
   printf "  Epsilon"
   printf "    G: Kin.E."
   printf "\n"
   for (i=1; i<=ncrit; i++) {
      if (type[i] != type_bnd  && type[i] != type_nuc) {
         epsilon = 100 * (1 - hess2[i] / hess3[i])
         printf "CP OTHER "
         printf "%-6s ", critpoint[type[i]]
         if (type[i] == type_nuc) {
            printf "%-4s(%02d) ", atname[neig1[i]], neig1[i]
            }
         else { printf "         " }
         printf "%-6s ", crittype[i]
         printf "%12.6f ", rho[i]
         if (type[i] != type_nuc) {
            printf "%12.6f ", del2[i]
            printf "%12.6f %8.2f ", hess1[i], epsilon
            printf "%12.6f ", kinG[i]
            }
         printf "\n"
         #MDC#if debug
            printf "coord: %12.6f %12.6f %12.6f", xcrit[i], ycrit[i], zcrit[i]
            printf "%-4s(%02d) ", atname[neig1[i]], neig1[i]
            printf "%-4s(%02d) ", atname[neig2[i]], neig2[i]
            printf "%9.3f %9.3f ", dist1[i], dist2[i]
            printf "\n"
         #MDC#endif debug
         }
      }
   printf "\n"

   printf "#Total number of c.p.: %4d\n", ncrit
   printf "#Types of c.p. (nbrc): %4d %4d %4d %4d\n", N_nuc, N_bnd, N_rng, N_min
   printf "\n"
   printf "#######################################################################\n"
   printf "\n"
   }
