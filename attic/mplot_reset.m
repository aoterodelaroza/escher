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

function mplot_reset ()
% function mplot_reset ()
%
% mplot_reset - restart the molecular plot database.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: December 2011

global ge 

ge = [];

ge.npts = 0;
ge.nstk = 0;
ge.nball = 0;
ge.bicolor = true;
ge.defstkrad = 0.1;
ge.defstkcolor = 121;

ge.jmlcolor = [ \
       255,255,255; 217,255,255; 204,128,255; # 001-003 H,  He, Li
       194,255,000; 255,181,181; 144,144,144; # 004-006 Be, B,  C
       048,080,248; 255,013,013; 144,224,080; # 007-009 N,  O,  F
       179,227,245; 171,092,242; 138,255,000; # 010-012 Ne, Na, Mg
       191,166,166; 240,200,160; 255,128,000; # 013-015 Al, Si, P
       255,255,048; 031,240,031; 128,209,227; # 016,018 S,  Cl, Ar
       143,064,212; 061,255,000; 230,230,230; # 019-021 K,  Ca, Sc
       191,194,199; 166,166,171; 138,153,199; # 022-024 Ti, V,  Cr
       156,122,199; 224,102,051; 240,144,160; # 025-027 Mn, Fe, Co
       080,208,080; 200,128,051; 125,128,176; # 028-030 Ni, Cu, Zn
       194,143,143; 102,143,143; 189,128,227; # 031-033 Ga, Ge, As
       255,161,000; 166,041,041; 092,184,209; # 034-036 Se, Br, Kr
       112,046,176; 000,255,000; 148,255,255; # 037-039 Rb, Sr, Y
       148,224,224; 115,194,201; 084,181,181; # 040-042 Zr, Nb, Mo
       059,158,158; 036,143,143; 010,125,140; # 043-045 Tc, Ru, Rh
       000,105,133; 192,192,192; 255,217,143; # 046-048 Pd, Ag, Cd
       166,117,115; 102,128,128; 158,099,181; # 049-051 In, Sn, Sb
       212,122,000; 148,000,148; 066,158,176; # 052-054 Te, I,  Xe
       087,023,143; 000,201,000; 112,212,255; # 055-057 Cs, Ba, La
       255,255,199; 217,255,199; 199,255,199; # 058-060 Ce, Pr, Nd
       163,255,199; 143,255,199; 097,255,199; # 061-063 Pm, Sm, Eu
       069,255,199; 048,255,199; 031,255,199; # 064-066 Gd, Tb, Dy
       000,255,156; 000,230,117; 000,212,082; # 067-069 Ho, Er, Tm
       000,191,056; 000,171,036; 077,194,255; # 070-072 Yb, Lu, Hf
       077,166,255; 033,148,214; 038,125,171; # 073-075 Ta, W,  Re
       038,102,150; 023,084,135; 208,208,224; # 076-078 Os, Ir, Pt
       255,209,035; 184,184,208; 166,084,077; # 079-081 Au, Hg, Tl
       087,089,097; 158,079,181; 171,092,000; # 082-084 Pb, Bi, Po
       117,079,069; 066,130,150; 066,000,102; # 085-087 At, Rn, Fr
       000,125,000; 112,171,250; 000,186,255; # 088-090 Ra, Ac, Th
       000,161,255; 000,143,255; 000,128,255; # 091-093 Pa, U,  Np
       000,107,255; 084,092,242; 120,092,227; # 094-096 Pu, Am, Cm
       138,079,227; 161,054,212; 179,031,212; # 097-099 Bk, Cf, Es
       179,031,186; 179,013,166; 189,013,135; # 100-102 Fm, Md, No
       199,000,102; 204,000,089; 209,000,079; # 103-105 Lr, Rf, Db
       217,000,069; 224,000,056; 230,000,046; # 106-108 Sg, Bh, Hs
       235,000,038; 160,000,066; 015,130,015; # 109-111 Mt, XX, YY
       020,090,255; 200,000,200; 255,180,070; # 112-114 ZZ, XY, XZ
       000,220,220; 230,010,010; 140,255,140; # 115-117 YX, YZ, ZX
       112,112,255; 127,050,000; 127,000,050; # 118-120 ZY, TT, WW
       ]';

endfunction
