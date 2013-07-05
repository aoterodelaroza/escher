#! /usr/bin/octave -q

addpath("../../../src");
cr = struct();
acell = 2;
cr.a = [acell acell acell];
cr.b = [90 90 90] * pi /180;
cr.nat = 8;
cr.ntyp = 2;
cr.x = [0   0.5 0.5
        0.5 0   0.5
        0.5 0.5 0
        0   0   0
        0.5 0.5 0.5
        0.5 0   0
        0   0.5 0
        0   0   0.5];
cr.typ = [1 1 1 1 2 2 2 2];
cr.qtyp = [1. -1.];
eps = 1d-15;
ene = cr_qewald(cr,eps);

printf("Madelung constant for NaCl, M = 1.7476\n")
printf(" M = -E_latt/4 for a cubic NaCl with a=2 and atomic charges = 1.\n")
printf("Total electrostatic energy : %.15f\n",-ene/4);

