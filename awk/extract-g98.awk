#! /usr/bin/awk -f

#
# extract-g98.awk - Analysis of the g98 output
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
   if (xx > dij || -xx > dij) { return (0) }
   else if (yy > dij || -yy > dij) { return (0) }
   else if (zz > dij || -zz > dij ) { return (0) }
   else if (xx*xx+yy*yy+zz*zz > dij*dij) { return (0) }
   else { return (1) }
   }

function PrintAngle(i1, i2, i3){
   x21 = atxc[i2,bng] - atxc[i1,bng]
   y21 = atyc[i2,bng] - atyc[i1,bng]
   z21 = atzc[i2,bng] - atzc[i1,bng]
   r21 = sqrt(x21*x21 + y21*y21 + z21*z21)
   x23 = atxc[i2,bng] - atxc[i3,bng]
   y23 = atyc[i2,bng] - atyc[i3,bng]
   z23 = atzc[i2,bng] - atzc[i3,bng]
   r23 = sqrt(x23*x23 + y23*y23 + z23*z23)
   pesc = (x21*x23 + y21*y23 + z21*z23)
   x123 = y21*z23 - y23*z21
   y123 = z21*x23 - z23*x21
   z123 = x21*y23 - x23*y21
   pvec = sqrt(x123^2 + y123^2 + z123^2)
   angle = atan2(pvec,pesc)
   angle = angle * Rad2Deg
   lbl0 = sprintf("%s-", atname[i1,bng])
   lbl0 = lbl0 sprintf("%s-", atname[i2,bng])
   lbl0 = lbl0 sprintf("%s", atname[i3,bng])
   lbl1 = sprintf("%s(%d)-", atname[i1,bng], i1)
   lbl1 = lbl1 sprintf("%s(%d)-", atname[i2,bng], i2)
   lbl1 = lbl1 sprintf("%s(%d)", atname[i3,bng], i3)
   printf "ANGLE %12.6f %-12s %s\n", angle, lbl0, lbl1
   }

function PrintDihed(i1, i2, i3, i4){
   # The vectorial products 12x23 and 23x34 will provide
   # the normal vectors of the two intersecting planes:
   x1 = atxc[i2,bng]-atxc[i1,bng]
   y1 = atyc[i2,bng]-atyc[i1,bng]
   z1 = atzc[i2,bng]-atzc[i1,bng]
   x2 = atxc[i3,bng]-atxc[i2,bng]
   y2 = atyc[i3,bng]-atyc[i2,bng]
   z2 = atzc[i3,bng]-atzc[i2,bng]
   ux1 = y1*z2-z1*y2
   uy1 = z1*x2-x1*z2
   uz1 = x1*y2-y1*x2
   x3 = atxc[i4,bng]-atxc[i3,bng]
   y3 = atyc[i4,bng]-atyc[i3,bng]
   z3 = atzc[i4,bng]-atzc[i3,bng]
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
   lbl0 = sprintf("%s-", atname[i1,bng])
   lbl0 = lbl0 sprintf("%s-", atname[i2,bng])
   lbl0 = lbl0 sprintf("%s-", atname[i3,bng])
   lbl0 = lbl0 sprintf("%s", atname[i4,bng])
   lbl1 = sprintf("%s(%d)-", atname[i1,bng], i1)
   lbl1 = lbl1 sprintf("%s(%d)-", atname[i2,bng], i2)
   lbl1 = lbl1 sprintf("%s(%d)-", atname[i3,bng], i3)
   lbl1 = lbl1 sprintf("%s(%d)", atname[i4,bng], i4)
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
      x21tmp = atxc[i2tmp,bng] - atxc[i1tmp,bng]
      y21tmp = atyc[i2tmp,bng] - atyc[i1tmp,bng]
      z21tmp = atzc[i2tmp,bng] - atzc[i1tmp,bng]
      x23tmp = atxc[i2tmp,bng] - atxc[i3tmp,bng]
      y23tmp = atyc[i2tmp,bng] - atyc[i3tmp,bng]
      z23tmp = atzc[i2tmp,bng] - atzc[i3tmp,bng]
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
   bondmax = 1.25
   atomicname[  1] = "H" ;   atomicnumber["H" ] =   1;   rcov[  1] =  53;
   atomicname[  2] = "He";   atomicnumber["HE"] =   2;   rcov[  2] =  70;
   atomicname[  3] = "Li";   atomicnumber["LI"] =   3;   rcov[  3] =  68;
   atomicname[  4] = "Be";   atomicnumber["BE"] =   4;   rcov[  4] =  35;
   atomicname[  5] = "B" ;   atomicnumber["B" ] =   5;   rcov[  5] =  83;
   atomicname[  6] = "C" ;   atomicnumber["C" ] =   6;   rcov[  6] =  68;
   atomicname[  7] = "N" ;   atomicnumber["N" ] =   7;   rcov[  7] =  68;
   atomicname[  8] = "O" ;   atomicnumber["O" ] =   8;   rcov[  8] =  68;
   atomicname[  9] = "F" ;   atomicnumber["F" ] =   9;   rcov[  9] =  64;
   atomicname[ 10] = "Ne";   atomicnumber["NE"] =  10;   rcov[ 10] =  70;
   atomicname[ 11] = "Na";   atomicnumber["NA"] =  11;   rcov[ 11] =  97;
   atomicname[ 12] = "Mg";   atomicnumber["MG"] =  12;   rcov[ 12] = 110;
   atomicname[ 13] = "Al";   atomicnumber["AL"] =  13;   rcov[ 13] = 135;
   atomicname[ 14] = "Si";   atomicnumber["SI"] =  14;   rcov[ 14] = 120;
   atomicname[ 15] = "P" ;   atomicnumber["P" ] =  15;   rcov[ 15] = 105;
   atomicname[ 16] = "S" ;   atomicnumber["S" ] =  16;   rcov[ 16] = 102;
   atomicname[ 17] = "Cl";   atomicnumber["CL"] =  17;   rcov[ 17] =  99;
   atomicname[ 18] = "Ar";   atomicnumber["AR"] =  18;   rcov[ 18] =  70;
   atomicname[ 19] = "K" ;   atomicnumber["K" ] =  19;   rcov[ 19] = 133;
   atomicname[ 20] = "Ca";   atomicnumber["CA"] =  20;   rcov[ 20] =  99;
   atomicname[ 21] = "Sc";   atomicnumber["SC"] =  21;   rcov[ 21] = 144;
   atomicname[ 22] = "Ti";   atomicnumber["TI"] =  22;   rcov[ 22] = 147;
   atomicname[ 23] = "V" ;   atomicnumber["V" ] =  23;   rcov[ 23] = 133;
   atomicname[ 24] = "Cr";   atomicnumber["CR"] =  24;   rcov[ 24] = 135;
   atomicname[ 25] = "Mn";   atomicnumber["MN"] =  25;   rcov[ 25] = 135;
   atomicname[ 26] = "Fe";   atomicnumber["FE"] =  26;   rcov[ 26] = 134;
   atomicname[ 27] = "Co";   atomicnumber["CO"] =  27;   rcov[ 27] = 133;
   atomicname[ 28] = "Ni";   atomicnumber["NI"] =  28;   rcov[ 28] = 150;
   atomicname[ 29] = "Cu";   atomicnumber["CU"] =  29;   rcov[ 29] = 152;
   atomicname[ 30] = "Zn";   atomicnumber["ZN"] =  30;   rcov[ 30] = 145;
   atomicname[ 31] = "Ga";   atomicnumber["GA"] =  31;   rcov[ 31] = 122;
   atomicname[ 32] = "Ge";   atomicnumber["GE"] =  32;   rcov[ 32] = 117;
   atomicname[ 33] = "As";   atomicnumber["AS"] =  33;   rcov[ 33] = 121;
   atomicname[ 34] = "Se";   atomicnumber["SE"] =  34;   rcov[ 34] = 122;
   atomicname[ 35] = "Br";   atomicnumber["BR"] =  35;   rcov[ 35] = 121;
   atomicname[ 36] = "Kr";   atomicnumber["KR"] =  36;   rcov[ 36] = 191;
   atomicname[ 37] = "Rb";   atomicnumber["RB"] =  37;   rcov[ 37] = 147;
   atomicname[ 38] = "Sr";   atomicnumber["SR"] =  38;   rcov[ 38] = 112;
   atomicname[ 39] = "Y" ;   atomicnumber["Y" ] =  39;   rcov[ 39] = 178;
   atomicname[ 40] = "Zr";   atomicnumber["ZR"] =  40;   rcov[ 40] = 157;
   atomicname[ 41] = "Nb";   atomicnumber["NB"] =  41;   rcov[ 41] = 148;
   atomicname[ 42] = "Mo";   atomicnumber["MO"] =  42;   rcov[ 42] = 147;
   atomicname[ 43] = "Tc";   atomicnumber["TC"] =  43;   rcov[ 43] = 135;
   atomicname[ 44] = "Ru";   atomicnumber["RU"] =  44;   rcov[ 44] = 140;
   atomicname[ 45] = "Rh";   atomicnumber["RH"] =  45;   rcov[ 45] = 145;
   atomicname[ 46] = "Pd";   atomicnumber["PD"] =  46;   rcov[ 46] = 150;
   atomicname[ 47] = "Ag";   atomicnumber["AG"] =  47;   rcov[ 47] = 159;
   atomicname[ 48] = "Cd";   atomicnumber["CD"] =  48;   rcov[ 48] = 169;
   atomicname[ 49] = "In";   atomicnumber["IN"] =  49;   rcov[ 49] = 163;
   atomicname[ 50] = "Sn";   atomicnumber["SN"] =  50;   rcov[ 50] = 146;
   atomicname[ 51] = "Sb";   atomicnumber["SB"] =  51;   rcov[ 51] = 146;
   atomicname[ 52] = "Te";   atomicnumber["TE"] =  52;   rcov[ 52] = 147;
   atomicname[ 53] = "I" ;   atomicnumber["I" ] =  53;   rcov[ 53] = 140;
   atomicname[ 54] = "Xe";   atomicnumber["XE"] =  54;   rcov[ 54] = 198;
   atomicname[ 55] = "Cs";   atomicnumber["CS"] =  55;   rcov[ 55] = 167;
   atomicname[ 56] = "Ba";   atomicnumber["BA"] =  56;   rcov[ 56] = 134;
   atomicname[ 57] = "La";   atomicnumber["LA"] =  57;   rcov[ 57] = 187;
   atomicname[ 58] = "Ce";   atomicnumber["CE"] =  58;   rcov[ 58] = 183;
   atomicname[ 59] = "Pr";   atomicnumber["PR"] =  59;   rcov[ 59] = 182;
   atomicname[ 60] = "Nd";   atomicnumber["ND"] =  60;   rcov[ 60] = 181;
   atomicname[ 61] = "Pm";   atomicnumber["PM"] =  61;   rcov[ 61] = 180;
   atomicname[ 62] = "Sm";   atomicnumber["SM"] =  62;   rcov[ 62] = 180;
   atomicname[ 63] = "Eu";   atomicnumber["EU"] =  63;   rcov[ 63] = 199;
   atomicname[ 64] = "Gd";   atomicnumber["GD"] =  64;   rcov[ 64] = 179;
   atomicname[ 65] = "Tb";   atomicnumber["TB"] =  65;   rcov[ 65] = 176;
   atomicname[ 66] = "Dy";   atomicnumber["DY"] =  66;   rcov[ 66] = 175;
   atomicname[ 67] = "Ho";   atomicnumber["HO"] =  67;   rcov[ 67] = 174;
   atomicname[ 68] = "Er";   atomicnumber["ER"] =  68;   rcov[ 68] = 173;
   atomicname[ 69] = "Tm";   atomicnumber["TM"] =  69;   rcov[ 69] = 172;
   atomicname[ 70] = "Yb";   atomicnumber["YB"] =  70;   rcov[ 70] = 194;
   atomicname[ 71] = "Lu";   atomicnumber["LU"] =  71;   rcov[ 71] = 172;
   atomicname[ 72] = "Hf";   atomicnumber["HF"] =  72;   rcov[ 72] = 157;
   atomicname[ 73] = "Ta";   atomicnumber["TA"] =  73;   rcov[ 73] = 143;
   atomicname[ 74] = "W" ;   atomicnumber["W" ] =  74;   rcov[ 74] = 137;
   atomicname[ 75] = "Re";   atomicnumber["RE"] =  75;   rcov[ 75] = 135;
   atomicname[ 76] = "Os";   atomicnumber["OS"] =  76;   rcov[ 76] = 137;
   atomicname[ 77] = "Ir";   atomicnumber["IR"] =  77;   rcov[ 77] = 132;
   atomicname[ 78] = "Pt";   atomicnumber["PT"] =  78;   rcov[ 78] = 150;
   atomicname[ 79] = "Au";   atomicnumber["AU"] =  79;   rcov[ 79] = 150;
   atomicname[ 80] = "Hg";   atomicnumber["HG"] =  80;   rcov[ 80] = 170;
   atomicname[ 81] = "Tl";   atomicnumber["TL"] =  81;   rcov[ 81] = 155;
   atomicname[ 82] = "Pb";   atomicnumber["PB"] =  82;   rcov[ 82] = 154;
   atomicname[ 83] = "Bi";   atomicnumber["BI"] =  83;   rcov[ 83] = 154;
   atomicname[ 84] = "Po";   atomicnumber["PO"] =  84;   rcov[ 84] = 168;
   atomicname[ 85] = "At";   atomicnumber["AT"] =  85;   rcov[ 85] = 170;
   atomicname[ 86] = "Rn";   atomicnumber["RN"] =  86;   rcov[ 86] = 240;
   atomicname[ 87] = "Fr";   atomicnumber["FR"] =  87;   rcov[ 87] = 200;
   atomicname[ 88] = "Ra";   atomicnumber["RA"] =  88;   rcov[ 88] = 190;
   atomicname[ 89] = "Ac";   atomicnumber["AC"] =  89;   rcov[ 89] = 188;
   atomicname[ 90] = "Th";   atomicnumber["TH"] =  90;   rcov[ 90] = 179;
   atomicname[ 91] = "Pa";   atomicnumber["PA"] =  91;   rcov[ 91] = 161;
   atomicname[ 92] = "U" ;   atomicnumber["U" ] =  92;   rcov[ 92] = 158;
   atomicname[ 93] = "Np";   atomicnumber["NP"] =  93;   rcov[ 93] = 155;
   atomicname[ 94] = "Pu";   atomicnumber["PU"] =  94;   rcov[ 94] = 153;
   atomicname[ 95] = "Am";   atomicnumber["AM"] =  95;   rcov[ 95] = 151;
   atomicname[ 96] = "Cm";   atomicnumber["CM"] =  96;   rcov[ 96] = 151;
   atomicname[ 97] = "Bk";   atomicnumber["BK"] =  97;   rcov[ 97] = 151;
   atomicname[ 98] = "Cf";   atomicnumber["CF"] =  98;   rcov[ 98] = 151;
   atomicname[ 99] = "Es";   atomicnumber["ES"] =  99;   rcov[ 99] = 151;
   atomicname[100] = "Fm";   atomicnumber["FM"] = 100;   rcov[100] = 151;
   atomicname[101] = "Md";   atomicnumber["MD"] = 101;   rcov[101] = 151;
   atomicname[102] = "No";   atomicnumber["NO"] = 102;   rcov[102] = 151;
   atomicname[103] = "Lw";   atomicnumber["LW"] = 103;   rcov[103] = 151;
   esMP2 = "no"
   esCI  = "no"
   }

/Default route:/ {
   getline; getline
   nroute
   route = ""
   while ($0 !~ /------*/) {
      nroute++
      route = route $0 "\n"
      getline
      if (nroute > 5) { break }
      }
   }

/l101.exe/ {
   getline
   getline
   runtitle = $0
   }

/Standard basis:/ { basisset = $0 }

/Standard orientation/ { ng++ }

/Standard orientation/,/Rotational constants/ {
   if (NF==6 && $1+0>0) {
      atn[ng]++
      i = atn[ng]
      atZ[i,ng]  = $2
      atxc[i,ng] = $4
      atyc[i,ng] = $5
      atzc[i,ng] = $6
      if ($2>0 && $2<104) { atname[i,ng] = atomicname[$2] }
      else                { atname[i,ng] = "XX" }
      }
   }

/Rotational constants/ {
   rotx[ng] = $4; roty[ng] = $5; rotz[ng] = $6; 
   }

/SCF Done:/ {
   energy[ng] = $5
   cycles[ng] = $8
   getline
   convg[ng] = $3
   virial[ng] = $6
   getline
   spin[ng] = $3
   }

/ EUMP2 / { esMP2 = "si"; energMP2[ng] = $NF }

/E\(CI\)=/       { esCI = "si"; energCI[ng] = $NF }
/E\(CI,SIZE\)=/  { esCI = "si"; energCIsize[ng] = $NF }

/Job cpu time:/ { $1=""; $2=""; $3=""; cpu = $0 }

{ lastline = $0 }

END {

   bng = ng

   if (esMP2 == "si") {
      beneg = energMP2[ng]
      for (i=1; i<=ng; i++) {
         if (energMP2[i] < beneg) { bng = i; beneg = energMP2[i] }
         }
      }
   else {
      beneg = energy[ng]
      for (i=1; i<=ng; i++) {
         if (energy[i] < beneg) { bng = i; beneg = energy[i] }
         }
      }

   printf "%s\n", runtitle
   printf "%s\n\n", route
   printf "Best energy (hartree)   : %18.10f\n", energy[bng]
   if (esMP2 == "si") {
      printf "MP2 energy  (hartree)   : %18.10f\n", energMP2[bng]
      }
   if (esCI == "si") {
      printf "CI      energy (hartree): %18.10f\n", energCI[bng]
      printf "CI+size energy (hartree): %18.10f\n", energCIsize[bng]
      }
   printf "Convergence             : %18.10f\n", convg[bng]
   printf "Virial (-V/T)           : %18.10f\n", virial[bng]
   printf "Spin S**2               : %18.10f\n", spin[bng]
   printf "Basis set line          : %s\n", basisset
   printf "Number of geometries    : %d\n", ng
   printf "Best geometries         : %d\n", bng
   printf "Analyzed file name      : %s\n", FILENAME
   printf "g98 cpu time            : %s\n", cpu
   printf "Last output line        : %s\n", lastline
   printf "Rotational constants GHz: %12.7f %12.7f %12.7f\n",
                                     rotx[bng], roty[bng], rotz[bng]
   printf "\n"
   if (atn[bng] <= 0) {
      print "---- Using initial coordinates ----------"
      print "---- Coordinates & distances in bohr ----"
      bng = 0
      }
   else {
      print "---- Using final coordinates ----------------"
      print "---- Coordinates & distances in Angstrom ----"
      }
   printf "    Atom  number       x              y              z\n"
   for (i=1; i<=atn[bng]; i++) {
      printf "%3d %-6s%4d  %15.10f%15.10f%15.10f\n", i, atname[i,bng],
             atZ[i,bng], atxc[i,bng], atyc[i,bng], atzc[i,bng]
      }
   printf "\n"

   ##for (i=1; i<=atn[bng]; i++) {
   ##   atrc[i] = sqrt(atxc[i,bng]^2 + atyc[i,bng]^2 + atzc[i,bng]^2)
   ##   }

   # Determine bonded atoms and print bond distances
   for (i=1; i<=atn[bng]; i++) {
      for (j=1; j<i; j++) {
         xx = atxc[i,bng] - atxc[j,bng]
         yy = atyc[i,bng] - atyc[j,bng]
         zz = atzc[i,bng] - atzc[j,bng]
         dcov = bondmax * (rcov[atZ[i,bng]] + rcov[atZ[j,bng]]) / 100
         bond[i,j] = isbonded(xx, yy, zz, dcov)
         bond[j,i] = bond[i,j]
         if (bond[i,j]) {
            i1 = i; i2 = j
            lbl0 = sprintf("%s-", atname[i1,bng])
            lbl0 = lbl0 sprintf("%s", atname[i2,bng])
            lbl1 = sprintf("%s(%d)-", atname[i1,bng], i1)
            lbl1 = lbl1 sprintf("%s(%d)", atname[i2,bng], i2)
            dij = sqrt(xx*xx + yy*yy + zz*zz)
            printf "BOND %12.6f %-9s %s\n", dij, lbl0, lbl1
            }
         }
      }
   printf "\n"

   # Print significative bond angles
   for (i=1; i<=atn[bng]; i++) {
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
   for (i=1; i<=atn[bng]; i++) {
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

   # Identify rings of connected atoms:
   # Algorithmic limitation: Instead of searching for rings of arbitrary
   # size we will only examine the occurrence of n-rings, where n: 3-7.
   # The algorithm used here calls for a generalization using recursivity.
   nring = 0
   nsupring = 0
   for (i1=1; i1<=atn[bng]; i1++) {
      for (i2=i1+1; i2<=atn[bng]; i2++) {
         if (bond[i1,i2]) {
            for (i3=i2+1; i3<=atn[bng]; i3++) {
               if (bond[i2,i3] && bond[i3,i1]) {
                  # A 3-ring
                  nring++
                  natr[nring] = 3
                  iatr[1,nring] = i1
                  iatr[2,nring] = i2
                  iatr[3,nring] = i3
                  }
               else if (bond[i2,i3]) {
                  for (i4=i3+1; i4<=atn[bng]; i4++) {
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
                        for (i5=i4+1; i5<=atn[bng]; i5++) {
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
                              for (i6=i4+1; i6<=atn[bng]; i6++) {
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
