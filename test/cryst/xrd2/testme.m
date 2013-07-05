#! /usr/bin/octave -q

cr = cr_read_espresso("pro_dl.scf.out");
cr.qtyp = [0 0 0 0];
[t ih pinfo] = cr_xrd(cr,"pro_dl",:,:,[5 30]);

cr = cr_read_espresso("pro_l.scf.out");
cr.qtyp = [0 0 0 0];
[t ih pinfo] = cr_xrd(cr,"pro_l",:,:,[5 30]);

