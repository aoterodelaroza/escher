#! /usr/bin/awk -f

#
# extract-xyz-promolden.awk - Analysis of the promolden output
# Get xyz files for the atoms in the output, plus independent xyz for
# critical points.
#
#-----------------------------------------------------------------------
# CopyRight (c): Victor Lua~na, Jun. 2012
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

function abs(x){ return (x<0 ? -x : x) }

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
   atomicnumber["H" ] =   1;    rcov[  1] =  53;   atomicname[  1] = "H" ;
   atomicnumber["HE"] =   2;    rcov[  2] =  70;   atomicname[  2] = "He";
   atomicnumber["LI"] =   3;    rcov[  3] =  68;   atomicname[  3] = "Li";
   atomicnumber["BE"] =   4;    rcov[  4] =  35;   atomicname[  4] = "Be";
   atomicnumber["B" ] =   5;    rcov[  5] =  83;   atomicname[  5] = "B" ;
   atomicnumber["C" ] =   6;    rcov[  6] =  68;   atomicname[  6] = "C" ;
   atomicnumber["N" ] =   7;    rcov[  7] =  68;   atomicname[  7] = "N" ;
   atomicnumber["O" ] =   8;    rcov[  8] =  68;   atomicname[  8] = "O" ;
   atomicnumber["F" ] =   9;    rcov[  9] =  64;   atomicname[  9] = "F" ;
   atomicnumber["NE"] =  10;    rcov[ 10] =  70;   atomicname[ 10] = "Ne";
   atomicnumber["NA"] =  11;    rcov[ 11] =  97;   atomicname[ 11] = "Na";
   atomicnumber["MG"] =  12;    rcov[ 12] = 110;   atomicname[ 12] = "Mg";
   atomicnumber["AL"] =  13;    rcov[ 13] = 135;   atomicname[ 13] = "Al";
   atomicnumber["SI"] =  14;    rcov[ 14] = 120;   atomicname[ 14] = "Si";
   atomicnumber["P" ] =  15;    rcov[ 15] = 105;   atomicname[ 15] = "P" ;
   atomicnumber["S" ] =  16;    rcov[ 16] = 102;   atomicname[ 16] = "S" ;
   atomicnumber["CL"] =  17;    rcov[ 17] =  99;   atomicname[ 17] = "Cl";
   atomicnumber["AR"] =  18;    rcov[ 18] =  70;   atomicname[ 18] = "Ar";
   atomicnumber["K" ] =  19;    rcov[ 19] = 133;   atomicname[ 19] = "K" ;
   atomicnumber["CA"] =  20;    rcov[ 20] =  99;   atomicname[ 20] = "Ca";
   atomicnumber["SC"] =  21;    rcov[ 21] = 144;   atomicname[ 21] = "Sc";
   atomicnumber["TI"] =  22;    rcov[ 22] = 147;   atomicname[ 22] = "Ti";
   atomicnumber["V" ] =  23;    rcov[ 23] = 133;   atomicname[ 23] = "V" ;
   atomicnumber["CR"] =  24;    rcov[ 24] = 135;   atomicname[ 24] = "Cr";
   atomicnumber["MN"] =  25;    rcov[ 25] = 135;   atomicname[ 25] = "Mn";
   atomicnumber["FE"] =  26;    rcov[ 26] = 134;   atomicname[ 26] = "Fe";
   atomicnumber["CO"] =  27;    rcov[ 27] = 133;   atomicname[ 27] = "Co";
   atomicnumber["NI"] =  28;    rcov[ 28] = 150;   atomicname[ 28] = "Ni";
   atomicnumber["CU"] =  29;    rcov[ 29] = 152;   atomicname[ 29] = "Cu";
   atomicnumber["ZN"] =  30;    rcov[ 30] = 145;   atomicname[ 30] = "Zn";
   atomicnumber["GA"] =  31;    rcov[ 31] = 122;   atomicname[ 31] = "Ga";
   atomicnumber["GE"] =  32;    rcov[ 32] = 117;   atomicname[ 32] = "Ge";
   atomicnumber["AS"] =  33;    rcov[ 33] = 121;   atomicname[ 33] = "As";
   atomicnumber["SE"] =  34;    rcov[ 34] = 122;   atomicname[ 34] = "Se";
   atomicnumber["BR"] =  35;    rcov[ 35] = 121;   atomicname[ 35] = "Br";
   atomicnumber["KR"] =  36;    rcov[ 36] = 191;   atomicname[ 36] = "Kr";
   atomicnumber["RB"] =  37;    rcov[ 37] = 147;   atomicname[ 37] = "Rb";
   atomicnumber["SR"] =  38;    rcov[ 38] = 112;   atomicname[ 38] = "Sr";
   atomicnumber["Y" ] =  39;    rcov[ 39] = 178;   atomicname[ 39] = "Y" ;
   atomicnumber["ZR"] =  40;    rcov[ 40] = 157;   atomicname[ 40] = "Zr";
   atomicnumber["NB"] =  41;    rcov[ 41] = 148;   atomicname[ 41] = "Nb";
   atomicnumber["MO"] =  42;    rcov[ 42] = 147;   atomicname[ 42] = "Mo";
   atomicnumber["TC"] =  43;    rcov[ 43] = 135;   atomicname[ 43] = "Tc";
   atomicnumber["RU"] =  44;    rcov[ 44] = 140;   atomicname[ 44] = "Ru";
   atomicnumber["RH"] =  45;    rcov[ 45] = 145;   atomicname[ 45] = "Rh";
   atomicnumber["PD"] =  46;    rcov[ 46] = 150;   atomicname[ 46] = "Pd";
   atomicnumber["AG"] =  47;    rcov[ 47] = 159;   atomicname[ 47] = "Ag";
   atomicnumber["CD"] =  48;    rcov[ 48] = 169;   atomicname[ 48] = "Cd";
   atomicnumber["IN"] =  49;    rcov[ 49] = 163;   atomicname[ 49] = "In";
   atomicnumber["SN"] =  50;    rcov[ 50] = 146;   atomicname[ 50] = "Sn";
   atomicnumber["SB"] =  51;    rcov[ 51] = 146;   atomicname[ 51] = "Sb";
   atomicnumber["TE"] =  52;    rcov[ 52] = 147;   atomicname[ 52] = "Te";
   atomicnumber["I" ] =  53;    rcov[ 53] = 140;   atomicname[ 53] = "I" ;
   atomicnumber["XE"] =  54;    rcov[ 54] = 198;   atomicname[ 54] = "Xe";
   atomicnumber["CS"] =  55;    rcov[ 55] = 167;   atomicname[ 55] = "Cs";
   atomicnumber["BA"] =  56;    rcov[ 56] = 134;   atomicname[ 56] = "Ba";
   atomicnumber["LA"] =  57;    rcov[ 57] = 187;   atomicname[ 57] = "La";
   atomicnumber["CE"] =  58;    rcov[ 58] = 183;   atomicname[ 58] = "Ce";
   atomicnumber["PR"] =  59;    rcov[ 59] = 182;   atomicname[ 59] = "Pr";
   atomicnumber["ND"] =  60;    rcov[ 60] = 181;   atomicname[ 60] = "Nd";
   atomicnumber["PM"] =  61;    rcov[ 61] = 180;   atomicname[ 61] = "Pm";
   atomicnumber["SM"] =  62;    rcov[ 62] = 180;   atomicname[ 62] = "Sm";
   atomicnumber["EU"] =  63;    rcov[ 63] = 199;   atomicname[ 63] = "Eu";
   atomicnumber["GD"] =  64;    rcov[ 64] = 179;   atomicname[ 64] = "Gd";
   atomicnumber["TB"] =  65;    rcov[ 65] = 176;   atomicname[ 65] = "Tb";
   atomicnumber["DY"] =  66;    rcov[ 66] = 175;   atomicname[ 66] = "Dy";
   atomicnumber["HO"] =  67;    rcov[ 67] = 174;   atomicname[ 67] = "Ho";
   atomicnumber["ER"] =  68;    rcov[ 68] = 173;   atomicname[ 68] = "Er";
   atomicnumber["TM"] =  69;    rcov[ 69] = 172;   atomicname[ 69] = "Tm";
   atomicnumber["YB"] =  70;    rcov[ 70] = 194;   atomicname[ 70] = "Yb";
   atomicnumber["LU"] =  71;    rcov[ 71] = 172;   atomicname[ 71] = "Lu";
   atomicnumber["HF"] =  72;    rcov[ 72] = 157;   atomicname[ 72] = "Hf";
   atomicnumber["TA"] =  73;    rcov[ 73] = 143;   atomicname[ 73] = "Ta";
   atomicnumber["W" ] =  74;    rcov[ 74] = 137;   atomicname[ 74] = "W" ;
   atomicnumber["RE"] =  75;    rcov[ 75] = 135;   atomicname[ 75] = "Re";
   atomicnumber["OS"] =  76;    rcov[ 76] = 137;   atomicname[ 76] = "Os";
   atomicnumber["IR"] =  77;    rcov[ 77] = 132;   atomicname[ 77] = "Ir";
   atomicnumber["PT"] =  78;    rcov[ 78] = 150;   atomicname[ 78] = "Pt";
   atomicnumber["AU"] =  79;    rcov[ 79] = 150;   atomicname[ 79] = "Au";
   atomicnumber["HG"] =  80;    rcov[ 80] = 170;   atomicname[ 80] = "Hg";
   atomicnumber["TL"] =  81;    rcov[ 81] = 155;   atomicname[ 81] = "Tl";
   atomicnumber["PB"] =  82;    rcov[ 82] = 154;   atomicname[ 82] = "Pb";
   atomicnumber["BI"] =  83;    rcov[ 83] = 154;   atomicname[ 83] = "Bi";
   atomicnumber["PO"] =  84;    rcov[ 84] = 168;   atomicname[ 84] = "Po";
   atomicnumber["AT"] =  85;    rcov[ 85] = 170;   atomicname[ 85] = "At";
   atomicnumber["RN"] =  86;    rcov[ 86] = 240;   atomicname[ 86] = "Rn";
   atomicnumber["FR"] =  87;    rcov[ 87] = 200;   atomicname[ 87] = "Fr";
   atomicnumber["RA"] =  88;    rcov[ 88] = 190;   atomicname[ 88] = "Ra";
   atomicnumber["AC"] =  89;    rcov[ 89] = 188;   atomicname[ 89] = "Ac";
   atomicnumber["TH"] =  90;    rcov[ 90] = 179;   atomicname[ 90] = "Th";
   atomicnumber["PA"] =  91;    rcov[ 91] = 161;   atomicname[ 91] = "Pa";
   atomicnumber["U" ] =  92;    rcov[ 92] = 158;   atomicname[ 92] = "U" ;
   atomicnumber["NP"] =  93;    rcov[ 93] = 155;   atomicname[ 93] = "Np";
   atomicnumber["PU"] =  94;    rcov[ 94] = 153;   atomicname[ 94] = "Pu";
   atomicnumber["AM"] =  95;    rcov[ 95] = 151;   atomicname[ 95] = "Am";
   atomicnumber["CM"] =  96;    rcov[ 96] = 151;   atomicname[ 96] = "Cm";
   atomicnumber["BK"] =  97;    rcov[ 97] = 151;   atomicname[ 97] = "Bk";
   atomicnumber["CF"] =  98;    rcov[ 98] = 151;   atomicname[ 98] = "Cf";
   atomicnumber["ES"] =  99;    rcov[ 99] = 151;   atomicname[ 99] = "Es";
   atomicnumber["FM"] = 100;    rcov[100] = 151;   atomicname[100] = "Fm";
   atomicnumber["MD"] = 101;    rcov[101] = 151;   atomicname[101] = "Md";
   atomicnumber["NO"] = 102;    rcov[102] = 151;   atomicname[102] = "No";
   atomicnumber["LW"] = 103;    rcov[103] = 151;   atomicname[103] = "Lw";

   bohr2angs = 0.52917720859
   angs2bohr = 1 / bohr2angs

   sym = ""
   title = "Unknown"
   mytitle = "Unknown"
   target = -100
   opt_step = -2
   doxyz = 0
   nxyz = 0
   dotopo = 0
   ntopo = 0
   ncp = 0
   typemax=1; typebnd=2; typerng=3; typemin=4; typeunk=5

   }

/File:/ { title = $NF; mytitle = title }

/Cartesian coordinates of centers:/ { target = FNR + 2; steps++; doxyz=1 }
FNR==target, /^ *$/ {
   if (NF>4) {
      nat = $1
      atZZ[nat] = atomicnumber[$2]
      atname[nat] = $2
      atxc[nat] = $3 * bohr2angs
      atyc[nat] = $4 * bohr2angs
      atzc[nat] = $5 * bohr2angs
      }
   else if (NF<=1) { target = -100 }
   }

doxyz==1 && FNR>target && /^ *$/ {
   file = sprintf("fxyz%03d.xyz", ++nxyz)
   printf " %d\n", nat > file
   printf "%s\n", mytitle >> file
   for (i=1; i<=nat; i++) {
      printf " %2s  %s %s %s\n", atname[i], atxc[i], atyc[i], atzc[i] >> file
      }
   close(file)
   doxyz = 0
   printf "%04d - %d atoms - %s\n", nxyz, nat, file
   }

/Analysis of Critical Point number/ { dotopo=1 }
###dotopo==1 && /CPoint at/
dotopo==1 && /CPoint at/ && /(3, -3)/ {
   tcp[++ncp]=typemax
   namecp[ncp]="max"
   namecp[ncp]="Cs"
   tx[ncp]=$(NF-2) * bohr2angs
   ty[ncp]=$(NF-1) * bohr2angs
   tz[ncp]=$(NF)   * bohr2angs
   for (i=1; i<=nat; i++) {
      tt = abs(atxc[i]-tx[ncp]) + abs(atyc[i]-ty[ncp]) + abs(atzc[i]-tz[ncp])
      if (tt <= 1e-5) { namecp[i] = atname[ncp] }
      }
   }
dotopo==1 && /CPoint at/ && /(3, -1)/ {
   tcp[++ncp]=typebnd
   namecp[ncp]="bnd"
   namecp[ncp]="Li"
   tx[ncp]=$(NF-2) * bohr2angs
   ty[ncp]=$(NF-1) * bohr2angs
   tz[ncp]=$(NF)   * bohr2angs
   }
dotopo==1 && /CPoint at/ && /(3,  1)/ {
   tcp[++ncp]=typerng
   namecp[ncp]="rng"
   namecp[ncp]="Be"
   tx[ncp]=$(NF-2) * bohr2angs
   ty[ncp]=$(NF-1) * bohr2angs
   tz[ncp]=$(NF)   * bohr2angs
   }
dotopo==1 && /CPoint at/ && /(3,  3)/ {
   tcp[++ncp]=typemin
   namecp[ncp]="min"
   namecp[ncp]="Xx"
   tx[ncp]=$(NF-2) * bohr2angs
   ty[ncp]=$(NF-1) * bohr2angs
   tz[ncp]=$(NF)   * bohr2angs
   }
dotopo==1 && /MGrad, Scal, Lapla/ {
   tscal[ncp]=$(NF-1)
   }
/ANALYSIS OF SYSTEM BONDS/ {
   dotopo=0

   # Critical points as an xyz file:
   mytitle = sprintf("Topology -- %s", mytitle)
   file = sprintf("txyz%03d.xyz", ++ntopo)
   printf " %d\n", ncp > file
   printf "Topology(MEP) %s\n", mytitle >> file
   for (i=1; i<=ncp; i++) {
      printf " %-2s  %s %s %s\n", namecp[i], tx[i], ty[i], tz[i] >> file
      }
   close(file)
   doxyz = 0
   printf "%04d - %d crit.p. - %s\n", ntopo, ncp, file

   # Critical points with additional data:
   # (not an xyz file, despite the initial name: changed to *.txt)
   ###file = sprintf("tmep%03d.xyz", ntopo)
   file = sprintf("tmep%03d.txt", ntopo)
   printf " %d\n", ncp > file
   printf "Topology(sMEP) %s\n", mytitle >> file
   for (i=1; i<=ncp; i++) {
      printf " %04d %-2s", i, namecp[i] >> file
      printf " %7.3f%7.3f%7.3f", tx[i], ty[i], tz[i] >> file
      printf " %10.6f %10.3f\n", tscal[i], tscal[i]*27.21138505 >> file
      }
   close(file)
   doxyz = 0
   printf "%04d - %d crit.p. - %s\n", ntopo, ncp, file

   # xyz file with only minima (lone pairs and similar):
   nminima = 0
   for (i=1;i<=ncp; i++) {
      if (tcp[i] == typemin) {++nminima}
   }
   file = sprintf("txyz-min.xyz")
   printf "%d\n", nminima > file
   printf "Topology(MEP) only MEP minima: %s\n", mytitle >> file
   for (i=1;i<=ncp; i++) {
       if (tcp[i] == typemin) {
          printf " Xx %s %s %s\n", tx[i], ty[i], tz[i] >> file
          }
       }
   close(file)
   doxyz = 0
   printf "%04d - %d minima - %s\n", ncp, nminima, file

   }
