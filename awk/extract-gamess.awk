#! /usr/bin/awk -f

#
# extract-gamess.awk - Analysis of the Gamess output
# Final geometry of an optimization and report of distances and angles
#
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

function isbonded(xx, yy, zz, dij){
   if      ( xx > dij) { return (0) }
   else if ( yy > dij) { return (0) }
   else if ( zz > dij) { return (0) }
   else if (-xx > dij) { return (0) }
   else if (-yy > dij) { return (0) }
   else if (-yy > dij) { return (0) }
   else if (xx*xx+yy*yy+zz*zz > dij*dij) { return (0) }
   else { return (1) }
   }

function PrintAngle(i1, i2, i3){
   x21 = atxc[i2,ng] - atxc[i1,ng]
   y21 = atyc[i2,ng] - atyc[i1,ng]
   z21 = atzc[i2,ng] - atzc[i1,ng]
   r21 = sqrt(x21*x21 + y21*y21 + z21*z21)
   x23 = atxc[i2,ng] - atxc[i3,ng]
   y23 = atyc[i2,ng] - atyc[i3,ng]
   z23 = atzc[i2,ng] - atzc[i3,ng]
   r23 = sqrt(x23*x23 + y23*y23 + z23*z23)
   pesc = (x21*x23 + y21*y23 + z21*z23)
   x123 = y21*z23 - y23*z21
   y123 = z21*x23 - z23*x21
   z123 = x21*y23 - x23*y21
   pvec = sqrt(x123^2 + y123^2 + z123^2)
   angle = atan2(pvec,pesc)
   angle = angle * Rad2Deg
   lbl0 = sprintf("%s-", atname[i1,ng])
   lbl0 = lbl0 sprintf("%s-", atname[i2,ng])
   lbl0 = lbl0 sprintf("%s", atname[i3,ng])
   lbl1 = sprintf("%s(%d)-", atname[i1,ng], i1)
   lbl1 = lbl1 sprintf("%s(%d)-", atname[i2,ng], i2)
   lbl1 = lbl1 sprintf("%s(%d)", atname[i3,ng], i3)
   printf "ANGLE %12.6f %-12s %s\n", angle, lbl0, lbl1
   }

function PrintDihed(i1, i2, i3, i4){
   # The vectorial products 12x23 and 23x34 will provide
   # the normal vectors of the two intersecting planes:
   x1 = atxc[i2,ng]-atxc[i1,ng]
   y1 = atyc[i2,ng]-atyc[i1,ng]
   z1 = atzc[i2,ng]-atzc[i1,ng]
   x2 = atxc[i3,ng]-atxc[i2,ng]
   y2 = atyc[i3,ng]-atyc[i2,ng]
   z2 = atzc[i3,ng]-atzc[i2,ng]
   ux1 = y1*z2-z1*y2
   uy1 = z1*x2-x1*z2
   uz1 = x1*y2-y1*x2
   x3 = atxc[i4,ng]-atxc[i3,ng]
   y3 = atyc[i4,ng]-atyc[i3,ng]
   z3 = atzc[i4,ng]-atzc[i3,ng]
   ux2 = z3*y2-y3*z2
   uy2 = x3*z2-z3*x2
   uz2 = y3*x2-x3*y2
   u1 = ux1*ux1+uy1*uy1+uz1*uz1
   u2 = ux2*ux2+uy2*uy2+uz2*uz2
   u = sqrt(u1*u2)
   if (u!=0.0) {
      cosa = (ux1*ux2+uy1*uy2+uz1*uz2) / u
      ux12 = uy1*uz2-uz1*uy2
      uy12 = uz1*ux2-ux1*uz2
      uz12 = ux1*uy2-uy1*ux2
      sina = sqrt(ux12^2 + uy12^2 + uz12^2) / u
      dihedr = atan2(sina,cosa) * Rad2Deg
      }
   else { dihedr = -720 }    # error!!!
   lbl0 = sprintf("%s-", atname[i1,ng])
   lbl0 = lbl0 sprintf("%s-", atname[i2,ng])
   lbl0 = lbl0 sprintf("%s-", atname[i3,ng])
   lbl0 = lbl0 sprintf("%s", atname[i4,ng])
   lbl1 = sprintf("%s(%d)-", atname[i1,ng], i1)
   lbl1 = lbl1 sprintf("%s(%d)-", atname[i2,ng], i2)
   lbl1 = lbl1 sprintf("%s(%d)-", atname[i3,ng], i3)
   lbl1 = lbl1 sprintf("%s(%d)", atname[i4,ng], i4)
   printf "DIHEDRAL %12.6f %-14s %s\n", dihedr, lbl0, lbl1
   }

function RingNormal(iring, normal){
   # Determines the average normal vector (components: xnor..znor)
   # for the ring number iring, and the average angular deviation from
   # planarity for the ring (angle angdev).
   nrtmp = natr[iring] - 2
   angdev = 0
   xn = 0
   yn = 0
   zn = 0
   # Decompose the ring in succesive triangles and average the normal
   # from all the triangles:
   for (itmp=1; itmp<=nrtmp; itmp++) {
      i1tmp = iatr[itmp,iring]
      i2tmp = iatr[itmp+1,iring]
      i3tmp = iatr[itmp+2,iring]
      x21tmp = atxc[i2tmp,ng] - atxc[i1tmp,ng]
      y21tmp = atyc[i2tmp,ng] - atyc[i1tmp,ng]
      z21tmp = atzc[i2tmp,ng] - atzc[i1tmp,ng]
      x23tmp = atxc[i2tmp,ng] - atxc[i3tmp,ng]
      y23tmp = atyc[i2tmp,ng] - atyc[i3tmp,ng]
      z23tmp = atzc[i2tmp,ng] - atzc[i3tmp,ng]
      x123tmp = y21tmp*z23tmp - y23tmp*z21tmp
      y123tmp = z21tmp*x23tmp - z23tmp*x21tmp
      z123tmp = x21tmp*y23tmp - x23tmp*y21tmp
      normtmp = sqrt(x123tmp^2 + y123tmp^2 + z123tmp^2)
      xntmp[itmp] = x123tmp / normtmp
      yntmp[itmp] = y123tmp / normtmp
      zntmp[itmp] = z123tmp / normtmp
      xn = xn + xntmp[itmp]
      yn = yn + yntmp[itmp]
      zn = zn + zntmp[itmp]
      }
   xn = xn / nrtmp
   yn = yn / nrtmp
   zn = zn / nrtmp
   normal[1] = xn
   normal[2] = yn
   normal[3] = zn
   for (itmp=1; itmp<=nrtmp; itmp++) {
      atmp = VecAngle(xn, yn, zn, xntmp[itmp], yntmp[itmp], zntmp[itmp])
      angdev = angdev + atmp
      }
   angdev = angdev / nrtmp
   return (angdev)
   }

function VecAngle(v1x, v1y, v1z, v2x, v2y, v2z){
   pesctmp = v1x*v2x + v1y*v2y + v1z*v2z
   xtmp = v1y*v2z - v1z*v2y
   ytmp = v1z*v2x - v1x*v2z
   ztmp = v1x*v2y - v1y*v2x
   pvectmp = sqrt(xtmp^2 + ytmp^2 + ztmp^2)
   atmpp = atan2(pvectmp,pesctmp)
   atmpp = atmpp * Rad2Deg
   return (atmpp)
   }


BEGIN{
   Rad2Deg = 45 / atan2(1,1)
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
   }

/RUN TITLE/ {
  getline
  getline
  gsub("^ *", "")
  gsub(" *$", "")
  runtitle = $0
  }

/STEP CPU/ { cpu = $10 }
/TOTAL WALL CLOCK/ { cpuutil = $10 }

/ATOM *ATOMIC *COORDINATES \(BOHR\)/,/^ *$/ {
   if (NF==5 && $2+0>0) {
      ng = 0
      atn[ng]++
      i = atn[ng]
      atname[i,ng] = $1
      atZ[i,ng]    = int($2)
      # Convert coordinates to Angstrom:
      atxc[i,ng]   = $3 * 0.5291771
      atyc[i,ng]   = $4 * 0.5291771
      atzc[i,ng]   = $5 * 0.5291771
      }
   }

/FINAL / && / ENERGY/ {
   ng++
   energy[ng] = $(NF-3)
   atn[ng] = 0
   }

/E\(MP2\)=/ {
   emp2 = $2
   }

/COORDINATES OF ALL ATOMS ARE/,/^$/ {
   if (NF==5 && $2+0>0) {
      atn[ng]++
      i = atn[ng]
      atname[i,ng] = $1
      atZ[i,ng]    = int($2)
      atxc[i,ng]   = $3
      atyc[i,ng]   = $4
      atzc[i,ng]   = $5
      }
   }

# TOTAL NUMBER OF SHELLS              =   27
# TOTAL NUMBER OF BASIS FUNCTIONS     =   64
# NUMBER OF ELECTRONS                 =   28
# CHARGE OF MOLECULE                  =    0
# STATE MULTIPLICITY                  =    1
# NUMBER OF OCCUPIED ORBITALS (ALPHA) =   14
# NUMBER OF OCCUPIED ORBITALS (BETA ) =   14
# TOTAL NUMBER OF ATOMS               =    8

END {

   printf "%s\n\n", runtitle
   printf "Final energy (hartree)  : %18.10f\n", energy[ng]
   if (emp2 < 0.0) {
      printf "MP2 energy   (hartree)  : %18.10f\n", emp2
      }
   if (ng > 1) {
      printf "Energy optimized        : %18.10f\n", energy[ng]-energy[1]
      printf "Last energy step        : %18.10f\n", energy[ng]-energy[ng-1]
      #for (i=2; i<=ng; i++) {
      #   printf "Step(%04d) energy change: %18.10f\n", i, energy[i]-energy[i-1]
      #   }
      }
   printf "Number of geometries    : %d\n", ng
   printf "Analyzed file name      : %s\n", FILENAME
   printf "gamess cpu time & usage : %d s (%s)\n", cpu, cpuutil
   printf "\n"
   if (atn[ng] <= 0) {
      print "---- Using initial coordinates ----------"
      print "---- Coordinates & distances in bohr ----"
      ng = 0
      }
   else {
      print "---- Using final coordinates ----------------"
      print "---- Coordinates & distances in Angstrom ----"
      }
   printf "    Atom  number       x              y              z\n"
   for (i=1; i<=atn[ng]; i++) {
      printf "%3d %-6s%4d  %15.10f%15.10f%15.10f\n", i, atname[i,ng],
             atZ[i,ng], atxc[i,ng], atyc[i,ng], atzc[i,ng]
      }
   printf "\n"

   ###for (i=1; i<=atn[ng]; i++) {
   ###   atrc[i] = sqrt(atxc[i,ng]^2 + atyc[i,ng]^2 + atzc[i,ng]^2)
   ###   }

   # Determine bonded atoms and print bond distances
   for (i=1; i<=atn[ng]; i++) {
      for (j=1; j<i; j++) {
         xx = atxc[i,ng] - atxc[j,ng]
         yy = atyc[i,ng] - atyc[j,ng]
         zz = atzc[i,ng] - atzc[j,ng]
         dcov = bondmax * (rcov[atZ[i,ng]] + rcov[atZ[j,ng]]) / 100
         bond[i,j] = isbonded(xx, yy, zz, dcov)
         bond[j,i] = bond[i,j]
         if (bond[i,j]) {
            i1 = i; i2 = j
            lbl0 = sprintf("%s-", atname[i1,ng])
            lbl0 = lbl0 sprintf("%s", atname[i2,ng])
            lbl1 = sprintf("%s(%d)-", atname[i1,ng], i1)
            lbl1 = lbl1 sprintf("%s(%d)", atname[i2,ng], i2)
            dij = sqrt(xx*xx + yy*yy + zz*zz)
            printf "BOND %12.6f %-9s %s\n", dij, lbl0, lbl1
            }
         }
      }
   printf "\n"

   # Print significative bond angles
   for (i=1; i<=atn[ng]; i++) {
      for (j=1; j<i; j++) {
         for (k=1; k<j; k++) {
            if (bond[i,j] && bond[j,k]) { PrintAngle(i,j,k) }
            if (bond[i,k] && bond[k,j]) { PrintAngle(i,k,j) }
            if (bond[j,i] && bond[i,k]) { PrintAngle(j,i,k) }
            }
         }
      }
   printf "\n"

   # Print significative dihedral angles
   for (i=1; i<=atn[ng]; i++) {
      for (j=1; j<i; j++) {
         for (k=1; k<j; k++) {
            for (l=1; l<k; l++) {
               # Possible connections:
               # kijl, lijk,   jikl, likj,   jilk, kilj
               # ijkl, ljki,   ijlk, kjli,   iklj, jkli
               if (bond[k,i] && bond[i,j] && bond[j,l]) { PrintDihed(k,i,j,l) }
               if (bond[l,i] && bond[i,j] && bond[j,k]) { PrintDihed(l,i,j,k) }
               if (bond[j,i] && bond[i,k] && bond[k,l]) { PrintDihed(j,i,k,l) }
               if (bond[l,i] && bond[i,k] && bond[k,j]) { PrintDihed(l,i,k,j) }
               if (bond[j,i] && bond[i,l] && bond[l,k]) { PrintDihed(j,i,l,k) }
               if (bond[k,i] && bond[i,l] && bond[l,j]) { PrintDihed(k,i,l,j) }
               if (bond[i,j] && bond[j,k] && bond[k,l]) { PrintDihed(i,j,k,l) }
               if (bond[l,j] && bond[j,k] && bond[k,i]) { PrintDihed(l,j,k,i) }
               if (bond[i,j] && bond[j,l] && bond[l,k]) { PrintDihed(i,j,l,k) }
               if (bond[k,j] && bond[j,l] && bond[l,i]) { PrintDihed(k,j,l,i) }
               if (bond[i,k] && bond[k,l] && bond[l,j]) { PrintDihed(i,k,l,j) }
               if (bond[j,k] && bond[k,l] && bond[l,i]) { PrintDihed(j,k,l,i) }
               }
            }
         }
      }

   exit(0)

   if (atn[ng] > 100) {
      print "Identification of ring structures skipped due to the big"
      print "number of atoms. The search could take ages!"
      exit(0)
      }

   # Identify rings of connected atoms:
   # Algorithmic limitation: Instead of searching for rings of arbitrary
   # size we will only examine the occurrence of n-rings, where n: 3--7.
   # The present algorithm can be improved and generalized using recursivity.
   nring = 0
   nsupring = 0
   for (i1=1; i1<=atn[ng]; i1++) {
      for (i2=i1+1; i2<=atn[ng]; i2++) {
         if (bond[i1,i2]) {
            for (i3=i2+1; i3<=atn[ng]; i3++) {
               if (bond[i2,i3] && bond[i3,i1]) {
                  # A 3-ring
                  nring++
                  natr[nring] = 3
                  iatr[1,nring] = i1
                  iatr[2,nring] = i2
                  iatr[3,nring] = i3
                  }
               else if (bond[i2,i3]) {
                  for (i4=i3+1; i4<=atn[ng]; i4++) {
                     if (bond[i3,i4] && bond[i4,i1]) {
                        # A 4-ring
                        nring++
                        natr[nring] = 4
                        iatr[1,nring] = i1
                        iatr[2,nring] = i2
                        iatr[3,nring] = i3
                        iatr[4,nring] = i4
                        }
                     else if (bond[i3,i4]) {
                        for (i5=i4+1; i5<=atn[ng]; i5++) {
                           if (bond[i4,i5] && bond[i5,i1]) {
                              # A 5-ring
                              nring++
                              natr[nring] = 5
                              iatr[1,nring] = i1
                              iatr[2,nring] = i2
                              iatr[3,nring] = i3
                              iatr[4,nring] = i4
                              iatr[5,nring] = i5
                              }
                           else if (bond[i4,i5]) {
                              for (i6=i4+1; i6<=atn[ng]; i6++) {
                                 if (bond[i5,i6] && bond[i6,i1]) {
                                    # A 6-ring
                                    nring++
                                    natr[nring] = 6
                                    iatr[1,nring] = i1
                                    iatr[2,nring] = i2
                                    iatr[3,nring] = i3
                                    iatr[4,nring] = i4
                                    iatr[5,nring] = i5
                                    iatr[6,nring] = i6
                                    }
                                 else if (bond[i5,i6]) {
                                    for (i7=i4+1; i7<=atn[ng]; i7++) {
                                       if (bond[i6,i7] && bond[i7,i1]) {
                                          # A 7-ring
                                          nring++
                                          natr[nring] = 7
                                          iatr[1,nring] = i1
                                          iatr[2,nring] = i2
                                          iatr[3,nring] = i3
                                          iatr[4,nring] = i4
                                          iatr[5,nring] = i5
                                          iatr[6,nring] = i6
                                          iatr[7,nring] = i7
                                          }
                                       else if (bond[i6,i7]) {
                                          # A possible ring with >7 atoms
                                          # Unsupported here
                                          nsupring++
                                          }
                                       } #for i7
                                    }
                                 } #for i6
                              }
                           } #for i5
                        }
                     } #for i4
                  }
               } #for i3
            }
         } #for i2
      } #for i1

   printf "\n"
   printf "Rings found in the structure : %6d\n", nring
   printf "Possible rings not determined: %6d\n", nsupring
   for (i=1; i<=nring; i++) {
      printf "  %3d) %d-ring formed by  : ", i, natr[i]
      for (j=1; j<=natr[i]; j++) {
         printf " %3d", iatr[j,i]
         }
      printf "\n"
      angdev = RingNormal(i, normal)
      printf "       nonplanarity angle: %12.6f\n", angdev
      printf "       normal vector     : "
      printf "%12.6f %12.6f %12.6f\n", normal[1], normal[2], normal[3]
      rnx[i] = normal[1]
      rny[i] = normal[2]
      rnz[i] = normal[3]
      }
   if (nring > 1) {
      printf "\n"
      printf "Dihedral angles between the rings:\n", nring
      for (i=1; i<=nring; i++) {
         for (j=i+1; j<=nring; j++) {
            printf "  Rings %3d and %3d", i, j
            ang = VecAngle(rnx[i], rny[i], rnz[i], rnx[j], rny[j], rnz[j])
            printf "  Angle %12.6f deg.\n", ang
            }
         }
      }

   }
