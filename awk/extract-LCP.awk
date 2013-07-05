#! /usr/bin/awk -f
#
# extract-LCP.awk - Extract Laplacian critical point (LCP) properties from
# the 'slightly verbose' output of extreme.
#
#-----------------------------------------------------------------------
# CopyRight (c): Victor Lua~na, May 20, 2000
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

function bubblesort(arr, isort, first, last,     idum, jdum, kdum){
      for (idum=first; idum<=last; idum++) {
         for (jdum=first; jdum<idum; jdum++) {
            if (arr[isort[jdum]] > arr[isort[idum]]) {
               kdum = isort[idum]
               isort[idum] = isort[jdum]
               isort[jdum] = kdum
               }
            }
         }
      }

function rainbow(x, xmin, xmax) {
   if (xmax != xmin) {
      rwx = 2 * (x-xmin) / (xmax-xmin) - 1
      }
   else {
      rwx = 1
      print "# ERROR en la escala arcoiris!"
      }
   rwred = 2 / (1 + exp((rwx-1)*(rwx-1)))
   rwgrn = 2 / (1 + exp(3 * rwx * rwx))
   rwblu = 2 / (1 + exp((rwx+1)*(rwx+1)))
   return( sprintf(" rgb %.4f %.4f %.4f", rwred, rwgrn, rwblu) )
   }

function arainbow(x, xmin, xmax) {
   if (xmax > xmin) {
      rwx1 = atan2(xmin,1e0)
      rwx2 = atan2(xmax,1e0)
      rwx = 2 * (atan2(x,1e0)-rwx1) / (rwx2-rwx1) - 1
      }
   else if (xmax < xmin) {
      rwx1 = atan2(xmax,1e0)
      rwx2 = atan2(xmin,1e0)
      rwx = 2 * (atan2(x,1e0)-rwx1) / (rwx2-rwx1) - 1
      }
   else {
      rwx = 1
      print "# ERROR en la escala arcoiris!"
      }
   #DBG# print "# x, xmin, xmax, rwx: ", x, xmin, xmax, rwx
   rwred = 2 / (1 + exp((rwx-1)*(rwx-1)))
   rwgrn = 2 / (1 + exp(3 * rwx * rwx))
   rwblu = 2 / (1 + exp((rwx+1)*(rwx+1)))
   return( sprintf(" rgb %.4f %.4f %.4f", rwred, rwgrn, rwblu) )
   }

function gray(x, xmin, xmax) {
   if (xmax != xmin) {
      grx = (xmax - x) / (xmax - xmin)
      }
   else {
      grx = 1
      print "# ERROR en la escala de grises!"
      }
   return( sprintf(" rgb %.4f %.4f %.4f", grx, grx, grx) )
   }

function abs(x) { return ( x >= 0 ? x : -x ) }


function AngleAlpha(i1, K2, i3){
   x21 = xcrit[K2] - atx[i1]
   y21 = ycrit[K2] - aty[i1]
   z21 = zcrit[K2] - atz[i1]
   r21 = sqrt(x21*x21 + y21*y21 + z21*z21)
   x23 = xcrit[K2] - atx[i3]
   y23 = ycrit[K2] - aty[i3]
   z23 = zcrit[K2] - atz[i3]
   r23 = sqrt(x23*x23 + y23*y23 + z23*z23)
   pesc = (x21*x23 + y21*y23 + z21*z23)
   x123 = y21*z23 - y23*z21
   y123 = z21*x23 - z23*x21
   z123 = x21*y23 - x23*y21
   pvec = sqrt(x123^2 + y123^2 + z123^2)
   angle = atan2(pvec,pesc)
   angle = 45 * angle / atan2(1,1)
   return (angle)
   }



BEGIN {
   MDIST = 6.0        # max distance for a non-bonded pair
   MinRho = 1e-6      # min value for the density
   bondmax = 1.15     # max elongation respect the sum of radii to have a bond
   bohr2a = 0.5291771
                       # Angstrom*100               Bohr
   atmZ["H" ] =   1;  rcov[  1] =  53;  rcore[  1] = 0.100;  col[  1] =   4
   atmZ["HE"] =   2;  rcov[  2] =  70;  rcore[  2] = 0.100;  col[  2] =   5
   atmZ["LI"] =   3;  rcov[  3] =  68;  rcore[  3] = 2.168;  col[  3] =  14
   atmZ["BE"] =   4;  rcov[  4] =  35;  rcore[  4] = 1.374;  col[  4] =  12
   atmZ["B" ] =   5;  rcov[  5] =  83;  rcore[  5] = 1.022;  col[  5] =  13
   atmZ["C" ] =   6;  rcov[  6] =  68;  rcore[  6] = 0.809;  col[  6] =  17
   atmZ["N" ] =   7;  rcov[  7] =  68;  rcore[  7] = 0.666;  col[  7] =  16
   atmZ["O" ] =   8;  rcov[  8] =  68;  rcore[  8] = 0.564;  col[  8] =   2
   atmZ["F" ] =   9;  rcov[  9] =  64;  rcore[  9] = 0.487;  col[  9] =  13
   atmZ["NE"] =  10;  rcov[ 10] =  70;  rcore[ 10] = 0.427;  col[ 10] =  12
   atmZ["NA"] =  11;  rcov[ 11] =  97;  rcore[ 11] = 3.194;  col[ 11] =   7
   atmZ["MG"] =  12;  rcov[ 12] = 110;  rcore[ 12] = 2.347;  col[ 12] =  15
   atmZ["AL"] =  13;  rcov[ 13] = 135;  rcore[ 13] = 1.890;  col[ 13] =   9
   atmZ["SI"] =  14;  rcov[ 14] = 120;  rcore[ 14] = 1.577;  col[ 14] =   6
   atmZ["P" ] =  15;  rcov[ 15] = 105;  rcore[ 15] = 1.354;  col[ 15] =   8
   atmZ["S" ] =  16;  rcov[ 16] = 102;  rcore[ 16] = 1.188;  col[ 16] =   3
   atmZ["CL"] =  17;  rcov[ 17] =  99;  rcore[ 17] = 1.057;  col[ 17] =  13
   atmZ["AR"] =  18;  rcov[ 18] =  70;  rcore[ 18] = 0.952;  col[ 18] =  12
   atmZ["K" ] =  19;  rcov[ 19] = 133;  rcore[ 19] = 0.863;  col[ 19] =  12
   atmZ["CA"] =  20;  rcov[ 20] =  99;  rcore[ 20] = 0.789;  col[ 20] =   9
   atmZ["SC"] =  21;  rcov[ 21] = 144;  rcore[ 21] = 0.733;  col[ 21] =  12
   atmZ["TI"] =  22;  rcov[ 22] = 147;  rcore[ 22] = 0.684;  col[ 22] =   9
   atmZ["V" ] =  23;  rcov[ 23] = 133;  rcore[ 23] = 0.642;  col[ 23] =  12
   atmZ["CR"] =  24;  rcov[ 24] = 135;  rcore[ 24] = 0.607;  col[ 24] =   9
   atmZ["MN"] =  25;  rcov[ 25] = 135;  rcore[ 25] = 0.571;  col[ 25] =   9
   atmZ["FE"] =  26;  rcov[ 26] = 134;  rcore[ 26] = 0.541;  col[ 26] =   8
   atmZ["CO"] =  27;  rcov[ 27] = 133;  rcore[ 27] = 0.514;  col[ 27] =  12
   atmZ["NI"] =  28;  rcov[ 28] = 150;  rcore[ 28] = 0.490;  col[ 28] =  10
   atmZ["CU"] =  29;  rcov[ 29] = 152;  rcore[ 29] = 0.469;  col[ 29] =  10
   atmZ["ZN"] =  30;  rcov[ 30] = 145;  rcore[ 30] = 0.447;  col[ 30] =  10
   atmZ["GA"] =  31;  rcov[ 31] = 122;  rcore[ 31] = 0.426;  col[ 31] =  12
   atmZ["GE"] =  32;  rcov[ 32] = 117;  rcore[ 32] = 0.407;  col[ 32] =  12
   atmZ["AS"] =  33;  rcov[ 33] = 121;  rcore[ 33] = 0.390;  col[ 33] =  12
   atmZ["SE"] =  34;  rcov[ 34] = 122;  rcore[ 34] = 0.374;  col[ 34] =  12
   atmZ["BR"] =  35;  rcov[ 35] = 121;  rcore[ 35] = 0.358;  col[ 35] =  10
   atmZ["KR"] =  36;  rcov[ 36] = 191;  rcore[ 36] = 1.417;  col[ 36] =  12
   atmZ["RB"] =  37;  rcov[ 37] = 147;  rcore[ 37] = 1.292;  col[ 37] =  12
   atmZ["SR"] =  38;  rcov[ 38] = 112;  rcore[ 38] = 1.195;  col[ 38] =  12
   atmZ["Y" ] =  39;  rcov[ 39] = 178;  rcore[ 39] = 1.117;  col[ 39] =  12
   atmZ["ZR"] =  40;  rcov[ 40] = 157;  rcore[ 40] = 1.050;  col[ 40] =  12
   atmZ["NB"] =  41;  rcov[ 41] = 148;  rcore[ 41] = 0.994;  col[ 41] =  12
   atmZ["MO"] =  42;  rcov[ 42] = 147;  rcore[ 42] = 0.942;  col[ 42] =  12
   atmZ["TC"] =  43;  rcov[ 43] = 135;  rcore[ 43] = 0.893;  col[ 43] =  12
   atmZ["RU"] =  44;  rcov[ 44] = 140;  rcore[ 44] = 0.853;  col[ 44] =  12
   atmZ["RH"] =  45;  rcov[ 45] = 145;  rcore[ 45] = 0.815;  col[ 45] =  12
   atmZ["PD"] =  46;  rcov[ 46] = 150;  rcore[ 46] = 0.781;  col[ 46] =  12
   atmZ["AG"] =  47;  rcov[ 47] = 159;  rcore[ 47] = 0.748;  col[ 47] =   9
   atmZ["CD"] =  48;  rcov[ 48] = 169;  rcore[ 48] = 0.717;  col[ 48] =  12
   atmZ["IN"] =  49;  rcov[ 49] = 163;  rcore[ 49] = 0.689;  col[ 49] =  12
   atmZ["SN"] =  50;  rcov[ 50] = 146;  rcore[ 50] = 0.662;  col[ 50] =  12
   atmZ["SB"] =  51;  rcov[ 51] = 146;  rcore[ 51] = 0.638;  col[ 51] =  12
   atmZ["TE"] =  52;  rcov[ 52] = 147;  rcore[ 52] = 0.615;  col[ 52] =  12
   atmZ["I" ] =  53;  rcov[ 53] = 140;  rcore[ 53] = 0.593;  col[ 53] =  11
   atmZ["XE"] =  54;  rcov[ 54] = 198;  rcore[ 54] = 0.573;  col[ 54] =  12
   atmZ["CS"] =  55;  rcov[ 55] = 167;  rcore[ 55] = 1.67 ;  col[ 55] =  12
   atmZ["BA"] =  56;  rcov[ 56] = 134;  rcore[ 56] = 1.34 ;  col[ 56] =   8
   atmZ["LA"] =  57;  rcov[ 57] = 187;  rcore[ 57] = 1.87 ;  col[ 57] =  12
   atmZ["CE"] =  58;  rcov[ 58] = 183;  rcore[ 58] = 1.83 ;  col[ 58] =  12
   atmZ["PR"] =  59;  rcov[ 59] = 182;  rcore[ 59] = 1.82 ;  col[ 59] =  12
   atmZ["ND"] =  60;  rcov[ 60] = 181;  rcore[ 60] = 1.81 ;  col[ 60] =  12
   atmZ["PM"] =  61;  rcov[ 61] = 180;  rcore[ 61] = 1.80 ;  col[ 61] =  12
   atmZ["SM"] =  62;  rcov[ 62] = 180;  rcore[ 62] = 1.80 ;  col[ 62] =  12
   atmZ["EU"] =  63;  rcov[ 63] = 199;  rcore[ 63] = 1.99 ;  col[ 63] =  12
   atmZ["GD"] =  64;  rcov[ 64] = 179;  rcore[ 64] = 1.79 ;  col[ 64] =  12
   atmZ["TB"] =  65;  rcov[ 65] = 176;  rcore[ 65] = 1.76 ;  col[ 65] =  12
   atmZ["DY"] =  66;  rcov[ 66] = 175;  rcore[ 66] = 1.75 ;  col[ 66] =  12
   atmZ["HO"] =  67;  rcov[ 67] = 174;  rcore[ 67] = 1.74 ;  col[ 67] =  12
   atmZ["ER"] =  68;  rcov[ 68] = 173;  rcore[ 68] = 1.73 ;  col[ 68] =  12
   atmZ["TM"] =  69;  rcov[ 69] = 172;  rcore[ 69] = 1.72 ;  col[ 69] =  12
   atmZ["YB"] =  70;  rcov[ 70] = 194;  rcore[ 70] = 1.94 ;  col[ 70] =  12
   atmZ["LU"] =  71;  rcov[ 71] = 172;  rcore[ 71] = 1.72 ;  col[ 71] =  12
   atmZ["HF"] =  72;  rcov[ 72] = 157;  rcore[ 72] = 1.57 ;  col[ 72] =  12
   atmZ["TA"] =  73;  rcov[ 73] = 143;  rcore[ 73] = 1.43 ;  col[ 73] =  12
   atmZ["W" ] =  74;  rcov[ 74] = 137;  rcore[ 74] = 1.37 ;  col[ 74] =  12
   atmZ["RE"] =  75;  rcov[ 75] = 135;  rcore[ 75] = 1.35 ;  col[ 75] =  12
   atmZ["OS"] =  76;  rcov[ 76] = 137;  rcore[ 76] = 1.37 ;  col[ 76] =  12
   atmZ["IR"] =  77;  rcov[ 77] = 132;  rcore[ 77] = 1.32 ;  col[ 77] =  12
   atmZ["PT"] =  78;  rcov[ 78] = 150;  rcore[ 78] = 1.50 ;  col[ 78] =  12
   atmZ["AU"] =  79;  rcov[ 79] = 150;  rcore[ 79] = 1.50 ;  col[ 79] =   6
   atmZ["HG"] =  80;  rcov[ 80] = 170;  rcore[ 80] = 1.70 ;  col[ 80] =  12
   atmZ["TL"] =  81;  rcov[ 81] = 155;  rcore[ 81] = 1.55 ;  col[ 81] =  12
   atmZ["PB"] =  82;  rcov[ 82] = 154;  rcore[ 82] = 1.54 ;  col[ 82] =  12
   atmZ["BI"] =  83;  rcov[ 83] = 154;  rcore[ 83] = 1.54 ;  col[ 83] =  12
   atmZ["PO"] =  84;  rcov[ 84] = 168;  rcore[ 84] = 1.68 ;  col[ 84] =  12
   atmZ["AT"] =  85;  rcov[ 85] = 170;  rcore[ 85] = 1.70 ;  col[ 85] =  12
   atmZ["RN"] =  86;  rcov[ 86] = 240;  rcore[ 86] = 2.40 ;  col[ 86] =  12
   atmZ["FR"] =  87;  rcov[ 87] = 200;  rcore[ 87] = 2.00 ;  col[ 87] =  12
   atmZ["RA"] =  88;  rcov[ 88] = 190;  rcore[ 88] = 1.90 ;  col[ 88] =  12
   atmZ["AC"] =  89;  rcov[ 89] = 188;  rcore[ 89] = 1.88 ;  col[ 89] =  12
   atmZ["TH"] =  90;  rcov[ 90] = 179;  rcore[ 90] = 1.79 ;  col[ 90] =  12
   atmZ["PA"] =  91;  rcov[ 91] = 161;  rcore[ 91] = 1.61 ;  col[ 91] =  12
   atmZ["U" ] =  92;  rcov[ 92] = 158;  rcore[ 92] = 1.58 ;  col[ 92] =  12
   atmZ["NP"] =  93;  rcov[ 93] = 155;  rcore[ 93] = 1.55 ;  col[ 93] =  12
   atmZ["PU"] =  94;  rcov[ 94] = 153;  rcore[ 94] = 1.53 ;  col[ 94] =  12
   atmZ["AM"] =  95;  rcov[ 95] = 151;  rcore[ 95] = 1.51 ;  col[ 95] =  12
   atmZ["CM"] =  96;  rcov[ 96] = 151;  rcore[ 96] = 1.51 ;  col[ 96] =  12
   atmZ["BK"] =  97;  rcov[ 97] = 151;  rcore[ 97] = 1.51 ;  col[ 97] =  12
   atmZ["CF"] =  98;  rcov[ 98] = 151;  rcore[ 98] = 1.51 ;  col[ 98] =  12
   atmZ["ES"] =  99;  rcov[ 99] = 151;  rcore[ 99] = 1.51 ;  col[ 99] =  12
   atmZ["FM"] = 100;  rcov[100] = 151;  rcore[100] = 1.51 ;  col[100] =  12
   atmZ["MD"] = 101;  rcov[101] = 151;  rcore[101] = 1.51 ;  col[101] =  12
   atmZ["NO"] = 102;  rcov[102] = 151;  rcore[102] = 1.51 ;  col[102] =  12
   atmZ["LW"] = 103;  rcov[103] = 151;  rcore[103] = 1.51 ;  col[103] =  12

   CPK[ 0] = "rgb 0.784 0.784 0.784"   # **0* LightGrey
   CPK[ 1] = "rgb 0.561 0.561 1.000"   # **1* SkyBlue
   CPK[ 2] = "rgb 0.941 0.000 0.000"   # **2* Red
   CPK[ 3] = "rgb 1.000 0.784 0.196"   # **3* Yellow
   CPK[ 4] = "rgb 1.000 1.000 1.000"   # **4* White
   CPK[ 5] = "rgb 1.000 0.753 0.796"   # **5* Pink
   CPK[ 6] = "rgb 0.855 0.647 0.125"   # **6* GoldenRed
   CPK[ 7] = "rgb 0.000 0.000 1.000"   # **7* Blue
   CPK[ 8] = "rgb 1.000 0.647 0.000"   # **8* Orange
   CPK[ 9] = "rgb 0.500 0.500 0.565"   # **9* DarkGrey
   CPK[10] = "rgb 0.647 0.165 0.165"   # *10* Brown
   CPK[11] = "rgb 0.627 0.078 0.567"   # *11* Purple
   CPK[12] = "rgb 1.000 0.078 0.567"   # *12* DeepPink
   CPK[13] = "rgb 0.000 0.900 0.000"   # *13* Green
   CPK[14] = "rgb 0.698 0.133 0.133"   # *14* FireBrick
   CPK[15] = "rgb 0.133 0.545 0.133"   # *15* MidGreen
   CPK[16] = "rgb 0.000 0.000 0.600"   # *16* DarkBlue
   CPK[17] = "rgb 0.600 0.600 0.700"   # *17* BluishMediumGrey

   eldens = 0
   eldel2rho = 0
   type_deg = -50
   type_nuc = -40
   type_nnm = -60
   type_sd1 = -30
   type_sd2 = -20
   type_min = -10
   type_atom = -99
   radius[type_atom] = 0.5
   radius[type_nnm] = 0.05
   radius[type_nuc] = 0.40
   radius[type_sd1] = 0.07
   radius[type_sd2] = 0.10
   radius[type_min] = 0.30
   radius[type_deg] = 0.03
   critpoint[type_deg] = "Deg?"
   critpoint[type_nuc] = "Nuc."
   critpoint[type_nnm] = "(3,-3)"
   critpoint[type_sd1] = "(3,-1)"
   critpoint[type_sd2] = "(3,+1)"
   critpoint[type_min] = "(3,+3)"
   ncp[type_deg] = 0
   ncp[type_nuc] = 0
   ncp[type_sd1] = 0
   ncp[type_sd2] = 0
   ncp[type_min] = 0
   ncp[type_nnm] = 0
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
      atiord[nat] = $2
      atx[nat] = $3
      aty[nat] = $4
      atz[nat] = $5
      atZnumber[nat] = atmZ[toupper(atname[nat])]
      }
   }

/^ ANALYZING RHO, THE ELECTRON DENSITY/ { eldens = 1 }
/^ ANALYZING DEL\*\*2(RHO), THE LAPLACIAN OF THE ELECTRON DENSITY/ {
   eldel2rho = 1 }

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
   if (rank < 3) { type[ncrit] = type_deg; ncp[type_deg]++ }
   else {
      if (signature == -3) { type[ncrit] = type_nuc; ncp[type_nuc]++ }
      else if (signature == -1) { type[ncrit] = type_sd1; ncp[type_sd1]++ }
      else if (signature == +1) { type[ncrit] = type_sd2; ncp[type_sd2]++ }
      else if (signature == +3) { type[ncrit] = type_min; ncp[type_min]++ }
      }
   crittype[ncrit] = sprintf("(%d,%d)", rank, signature)
   #
   # Sort the neighbors and guess to which atoms is this LCP connected:
   #
   bubblesort(dist, iord, 1, nneig,   idum, jdum, kdum) 
   neig1[ncrit] = iord[1]
   neig2[ncrit] = iord[2]
   neig3[ncrit] = iord[3]
   neig4[ncrit] = iord[4]
   dist1[ncrit] = dist[iord[1]]
   dist2[ncrit] = dist[iord[2]]
   dist3[ncrit] = dist[iord[3]]
   dist4[ncrit] = dist[iord[4]]
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
      if (dist1[ncrit] > 1e-2) {
         ncp[type_nuc]--
         ncp[type_nnm]++
         type[ncrit] = type_nnm
         }
      }
   }

/^ Rho\(r\)/             { rho[ncrit] = $2 }
/^ Del\*\*2\*\(Rho\(r\)\)/ { del2[ncrit] = $2 }
/^ DEL\*\*2\(DEL\*\*2\(Rho\(r\)\)\)/ { del2del2[ncrit] = $2 }
/^ G\(r\)/               { kinG[ncrit] = $2 }
/^ V\(r\)/               {
   #
   # Check if this is a valid critical point and remove it otherwise
   # Things to check:
   #                  is rho >= MinRho ?
   #                  is dist1[] <= MDIST ?
   #                  is dist1[] >= rcore[] ?
   #                  if this point repeated ?
   #
   remove = 0
   if (rho[ncrit] < MinRho) remove = 1
   if (dist1[ncrit] > MDIST) remove = 1
   iat = atmZ[toupper(atname[neig1[ncrit]])]
   rcor = rcore[iat]
   if (dist1[ncrit] < rcor) remove = 1
   i = 1
   while (remove == 0  &&  i < ncrit) {
      if ( abs(xcrit[ncrit] - xcrit[i]) < 1e-6  &&
           abs(ycrit[ncrit] - ycrit[i]) < 1e-6  &&
           abs(zcrit[ncrit] - zcrit[i]) < 1e-6  &&
           type[ncrit] == type[i] ) { remove = 1; N_repeated++ }
      else { i++ }
      }
   if (remove == 1) {
      N_removed++
      ncp[type[ncrit]]--
      ncrit--
      }
   }

END {
   #
   # Translate the origin to the center of mass ?
   #
   Xcm = 0
   Ycm = 0
   Zcm = 0
   for (i=1; i<=nat; i++) {
       Xcm += atx[i]
       Ycm += aty[i]
       Zcm += atz[i]
       }
   Xcm /= nat
   Ycm /= nat
   Zcm /= nat
   #
   # Start report. Print atomic positions:
   #
   printf "Datos extraidos del fichero: %s\n\n", FILENAME
   printf "%s\n\n", runtitle
   printf "Geometry center: %11.6f %11.6f %11.6f\n", Xcm, Ycm, Zcm
   printf "Atoms and positions\n"
   for (i=1; i<=nat; i++) {
      printf "%3d ", i
      printf "%-6s", atname[i]
      printf "%11.6f", atx[i]
      printf "%11.6f", aty[i]
      printf "%11.6f", atz[i]
      printf "\n"
      }
   printf "\n"

   printf "Geometrical distances between atoms\n"
   # Determine bonded atoms and print bond distances
   for (i=1; i<=nat; i++) {
      for (j=1; j<i; j++) {
         xx = atx[i] - atx[j]
         yy = aty[i] - aty[j]
         zz = atz[i] - atz[j]
         dij = sqrt(xx*xx + yy*yy + zz*zz)
         dcov = (rcov[atZnumber[i]] + rcov[atZnumber[j]]) / 100
         dcov = dcov / bohr2a
         bond[i,j] = (dij <= bondmax * dcov)
         bond[j,i] = bond[i,j]
         if (bond[i,j]) {
            i1 = i; i2 = j
            lbl0 = sprintf("%s-", atname[i1])
            lbl0 = lbl0 sprintf("%s", atname[i2])
            lbl1 = sprintf("%s(%d)-", atname[i1], i1)
            lbl1 = lbl1 sprintf("%s(%d)", atname[i2], i2)
            printf "BOND %12.6f %-9s %s\n", dij, lbl0, lbl1
            #
            # add to the list of neighbors
            neig_n[i]++
            neig_n[j]++
            neig_id[i,neig_n[i]] = j
            neig_id[j,neig_n[j]] = i
            neig_dis[i,neig_n[i]] = dij
            neig_dis[j,neig_n[j]] = dij
            }
         }
      }
   printf "\n"

   #
   # Sort critical points. First by type and within each type by distance
   # to the first neighbor:
   #
   inum = 0
   iold = 1
   for (i=1; i<=ncrit; i++) {
      if (type[i] == type_min) { inum++; iord[inum] = i }
      }
   bubblesort(dist1, iord, iold, inum,   idum, jdum, kdum) 
   iold = inum + 1
   for (i=1; i<=ncrit; i++) {
      if (type[i] == type_sd2) { inum++; iord[inum] = i }
      }
   bubblesort(dist1, iord, iold, inum,   idum, jdum, kdum) 
   iold = inum + 1
   for (i=1; i<=ncrit; i++) {
      if (type[i] == type_sd1) { inum++; iord[inum] = i }
      }
   bubblesort(dist1, iord, iold, inum,   idum, jdum, kdum) 
   iold = inum + 1
   for (i=1; i<=ncrit; i++) {
      if (type[i] == type_nnm) { inum++; iord[inum] = i }
      }
   bubblesort(dist1, iord, iold, inum,   idum, jdum, kdum) 
   iold = inum + 1
   for (i=1; i<=ncrit; i++) {
      if (type[i] == type_nuc) { inum++; iord[inum] = i }
      }
   bubblesort(dist1, iord, iold, inum,   idum, jdum, kdum) 
   iold = inum + 1
   for (i=1; i<=ncrit; i++) {
      if (type[i] == type_deg) { inum++; iord[inum] = i }
      }
   bubblesort(dist1, iord, iold, inum,   idum, jdum, kdum) 

   #
   # Report critical point position and properties:
   #
   printf "Other critical points\n"
   printf "    Num "
   printf "Type   "
   printf " Critical Point position (x,y,z)"
   printf "   E. density"
   printf "    Laplacian"
   printf "    G: Kin.E."
#   printf "          Hessian matrix eigenvalues     "
   printf "   Atom   Rad1 "
   printf "   Atom   Rad2 "
   printf "   Atom   Rad3 "
   printf "   Atom   Rad4 "
   printf "\n\n"
   typold = -9999
   for (i1=1; i1<=ncrit; i1++) {
      i = iord[i1]
      if (typold != type[i]) {for (k=1;k<=78;k++) {printf "-"}; printf "\n"}
      typold = type[i]
      printf "%3d %3d ", i1, i
      printf "%-6s", critpoint[type[i]]
      printf "%11.6f%11.6f%11.6f ", xcrit[i], ycrit[i], zcrit[i]
      printf "%12.6f ", rho[i]
      printf "%12.6f %12.6f ", del2[i], kinG[i]
#      printf "%12.6f %12.6f %12.6f ", hess1[i], hess2[i], hess3[i]
      #
      # report nearest neighbor (n1) and distance to it, plus analyze the
      # angles:
      #           alpha[i] :  n1 - LCP - nn[i]
      # where nn[i] is any of the atoms bonded to n1.
      #
      LonePair = 1
      n1 = neig1[i]
      printf "%4s(%02d)%7.4f", atname[neig1[i]], neig1[i], dist1[i]
      for (k = 1; k <= neig_n[n1]; k++) {
         n2 = neig_id[n1,k]
         alpha = AngleAlpha(n1, i, n2)
         ###DB##G#printf " DBG %4s(%02d)%7.4f", atname[n2], n2, alpha
         if (abs(alpha-180) < 15) {
            LonePair = 0
            printf " BOND-LCP %4s(%02d)%7.4f", atname[n2], n2, alpha
            }
         }
      if (LonePair == 1) { printf " LonePair" }
      #printf "%4s(%02d)%7.4f", atname[neig2[i]], neig2[i], dist2[i]
      #printf "%4s(%02d)%7.4f", atname[neig3[i]], neig3[i], dist3[i]
      #printf "%4s(%02d)%7.4f", atname[neig4[i]], neig4[i], dist4[i]
      printf "\n"
      if (i1%5 == 0) printf "\n"
      }
   printf "\n"

   printf "Total number of c.p.: %4d\n", ncrit
   printf "c.p. (Max,Nnm,Saddle,saddle,min,Deg): %4d %4d %4d %4d %4d %4d\n",
          ncp[type_nuc], ncp[type_nnm], ncp[type_sd1],
          ncp[type_sd2], ncp[type_min], ncp[type_deg]
   printf "Removed c.p. %d of which %d were repeated\n",
          N_removed, N_repeated

   #
   # Min & max rho values for each type of c.p.
   #
   for (i in ncp) {
      lapmin[i] = 1e60;  lapmax[i] = -1e60
      rhomin[i] = 1e60;  rhomax[i] = -1e60
      }
   for (i=1; i<=ncrit; i++) {
      iat = atmZ[toupper(atname[neig1[i]])]
      rcor = rcore[iat]
      if (dist1[i] >= rcor) {
         if (rho[i] < rhomin[type[i]]) rhomin[type[i]] = rho[i]
         if (rho[i] > rhomax[type[i]]) rhomax[type[i]] = rho[i]
         if (del2[i] < lapmin[type[i]]) lapmin[type[i]] = del2[i]
         if (del2[i] > lapmax[type[i]]) lapmax[type[i]] = del2[i]
         }
      }
   for (i in ncp) {
      if (ncp[i] > 0) {
         printf "Min & max density   for %s c.p.: %12.5e %12.5e\n",
                critpoint[i], rhomin[i], rhomax[i]
         printf "Min & max laplacian for %s c.p.: %12.5e %12.5e\n",
                critpoint[i], lapmin[i], lapmax[i]
         }
      }

   #
   # Output coordinate file for the MOLDEN visualization code
   #
   file = "coord-nnm.molden"
   print "[Molden Format]"  > file
   print "[Atoms] AU" > file
#   print nat+ncp[type_min] > file
#   print "Atomos y puntos (3,+3) en -del2rho" > file
   k = 0
   for (i=1; i<=nat; i++) {
      k++
      printf " %-4s %4d %6d %15.9f %15.9f %15.9f\n",
             atname[i], k, atmZ[toupper(atname[i])],
             atx[i], aty[i], atz[i]   > file
      }
   for (i=1; i<=ncrit; i++) {
      iat = atmZ[toupper(atname[neig1[i]])]
      rcor = rcore[iat]
      if (type[i] == type_min  &&  dist1[i] >= rcor) {
         k++
         printf " %-4s %4d %6d %15.9f %15.9f %15.9f\n",
                "XX", k, 99, xcrit[i], ycrit[i], zcrit[i]   > file
         }
      if (type[i] == type_sd2  &&  dist1[i] >= rcor) {
         k++
         printf " %-4s %4d %6d %15.9f %15.9f %15.9f\n",
                "YY", k, 100, xcrit[i], ycrit[i], zcrit[i]  > file
         }
      }
   close(file)

   #
   # Output file for the TESSEL code
   #
   for (i=1; i<=nat; i++) {
      atx[i] -= Xcm
      aty[i] -= Ycm
      atz[i] -= Zcm
      }
   for (i=1; i<=ncrit; i++) {
      xcrit[i] -= Xcm
      ycrit[i] -= Ycm
      zcrit[i] -= Zcm
      }
   file = "grafo.tess"
   print "#Atoms and Laplacian (3,+3) crit. points."            > file
   print "#set camangle 75 30 35"                               > file
   print "set camangle 25 -5 40"                                > file
   print "set background background {color rgb <0.6,0.9,0.8>}"  > file
   print "set use_planes .false."                               > file
   print "set zoom 2"                                           > file
   print "molecule"                                             > file
   print "   crystal"                                           > file
   print "      spg p 1"                                        > file
   print "      cell 1 1 1 90 90 90"                            > file
   for (i=1; i<=nat; i++) {
      printf "      NEQCLUS 1 %-4s\n", atname[i]                > file
      printf "      %15.9f %15.9f %15.9f\n",
                                       atx[i], aty[i], atz[i]   > file
      atomos[atname[i]] = atmZ[toupper(atname[i])]
      }
   for (i=1; i<=ncrit; i++) {
      iat = atmZ[toupper(atname[neig1[i]])]
      rcor = rcore[iat]
      if (type[i] == type_min  &&  dist1[i] >= rcor) {
         printf "      NEQCLUS 1 XX%03d\n", i                   > file
         printf "      %15.9f %15.9f %15.9f\n",
                                 xcrit[i], ycrit[i], zcrit[i]   > file
         }
      if (type[i] == type_sd2  &&  dist1[i] >= rcor) {
         printf "#      NEQCLUS 1 YY%03d\n", i                  > file
         printf "#      %15.9f %15.9f %15.9f\n",
                                 xcrit[i], ycrit[i], zcrit[i]   > file
         }
      }
   print "      crystalbox -100 -100 -100 100 100 100"          > file
   print "      clippingbox -100 -100 -100 100 100 100"         > file
   print "   endcrystal"                                        > file
   for (i in atomos) {
      iZ = atomos[i]
      printf "   BALL %-4s radius %7.3f %s\n",
                    i, radius[type_atom], CPK[col[iZ]]  > file
      }
   for (i in atomos) {
      for (j in atomos) {
         ###print i, j, atomos[i], atomos[j]
         if (atomos[i] > atomos[j]) {
            printf "   STICK %-4s %-4s radius 0.05\n", i, j     > file
            }
         }
      }
   for (i=1; i<=ncrit; i++) {
      iat = atmZ[toupper(atname[neig1[i]])]
      rcor = rcore[iat]
      if ( type[i] == type_min && dist1[i] >= rcor ) {
         it = type[i]
         printf "   BALL XX%03d radius %7.3f %s\n",
                    i, radius[type[i]],
                    arainbow(-del2[i],-lapmax[it],-lapmin[it])  > file
         }
      }
   for (i=1; i<=ncrit; i++) {
      iat = atmZ[toupper(atname[neig1[i]])]
      rcor = rcore[iat]
      if ( type[i] == type_sd2 && dist1[i] >= rcor ) {
         it = type[i]
         printf "#   BALL YY%03d radius %7.3f %s\n",
                    i, radius[type[i]],
                    arainbow(-del2[i],-lapmax[it],-lapmin[it])  > file
         }
      }
   for (i=1; i<=ncrit; i++) {
      iat = atmZ[toupper(atname[neig1[i]])]
      rcor = rcore[iat]
      if (type[i] == type_min  &&  dist1[i] >= rcor) {
         it = type[i]
         printf "   STICK %s XX%03d DIST %.4f %.4f RADIUS %7.3f %s\n",
                    atname[neig1[i]], i, dist1[i]*0.9, dist1[i]*1.1,
                    radius[it]*0.8,
                    rainbow(rho[i],rhomin[it],rhomax[it])       > file
         }
      }
   for (i=1; i<=ncrit; i++) {
      iat = atmZ[toupper(atname[neig1[i]])]
      rcor = rcore[iat]
      if (type[i] == type_sd2  &&  dist1[i] >= rcor) {
         it = type[i]
         printf "#   STICK %s YY%03d DIST %.4f %.4f RADIUS %7.3f %s\n",
                    atname[neig1[i]], i, dist1[i]*0.9, dist1[i]*1.1,
                    radius[it]*0.8,
                    rainbow(rho[i],rhomin[it],rhomax[it])       > file
         printf "#   STICK %s YY%03d DIST %.4f %.4f RADIUS %7.3f %s\n",
                    atname[neig2[i]], i, dist2[i]*0.9, dist2[i]*1.1,
                    radius[it]*0.8,
                    rainbow(rho[i],rhomin[it],rhomax[it])       > file
         }
      }
   print "  povray grafo.pov"                                   > file
   print "#  vrml grafo.vrml"                                   > file
   print "  off grafo.off"                                      > file
   print "endmolecule"                                          > file
   print "#run povray grafo.pov +W100 +H100"                    > file
   print "#run povray grafo.pov +W400 +H400 +A0.2 +R3"          > file
   print "#run pov3 grafo.pov +W100 +H100"                      > file
   print "run pov3 grafo.pov +W600 +H600 +A"                    > file
   print "run cjpeg -Q 80 grafo.tga > grafo.jpg"                > file
   print "run rm grafo.tga"                                     > file
   print "reset"                                                > file
   print "end"                                                  > file

   }
