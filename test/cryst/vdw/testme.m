#! /usr/bin/octave -q

addpath("../../../src");

## fcc Ar
## brute-force summation gives -1.484706724122E-02 Hy
cr = struct();
acell = 10;
cr.a = [acell acell acell];
cr.b = [90 90 90] * pi /180;
cr.nat = 4;
cr.ntyp = 1;
cr.x = [0 0.5 0.5
        0.5 0 0.5
        0.5 0.5 0
        0 0 0];
cr.typ = [1 1 1 1];
cr.c6typ = zeros(1);
cr.c6typ(1,1) = 64.2;
eps = 1d-15;
ene = cr_vdwewald(cr,eps);
printf("Edisp (Ar fcc, direct sum) : %.10e\n",-1.484706724122E-02);
printf("Edisp (Ar fcc, Ewald)      : %.10e\n",ene);

## graphite
cr.a = [4.641167  4.641167  12.000000];
cr.b = [90 90 120] * pi / 180;
cr.nat = 4;
cr.ntyp = 1;
cr.x = [0 0 1/4
        0 0 3/4
        2/3 1/3 1/4
        1/3 2/3 3/4];
cr.typ = [1 1 1 1];
cr.c6typ = [46.6];
cr.rvdwtyp = [6.5];
eps = 1d-15;
ene = cr_vdwewald(cr,eps);
printf("Edisp (graphite, no damping) : %.10e\n",ene);
ene = cr_vdwewald(cr,eps,1,[20 0.94]);
printf("Edisp (graphite, fermi damping) : %.10e\n",ene);
ene = cr_vdwewald(cr,eps,2,[0.8 1.0]);
printf("Edisp (graphite, bj damping) : %.10e\n",ene);

## cdcl2, high p
cr.a = [8.486018 8.486018 8.486017];
cr.b = [38.558081 38.558081 38.558081] * pi / 180;
cr.nat = 3;
cr.ntyp = 2;
cr.x = [ 0.0000000D+00  0.0000000D+00  0.0000000D+00
         0.2286240D+00  0.2286240D+00  0.2286240D+00
         0.7713760D+00  0.7713760D+00  0.7713760D+00];
cr.typ = [1 2 2];
cr.c6typ = [ 466.000000 194.238745
	     194.238745 94.600000];
rvdw = [4.043436 3.709223];
eps = 1d-15;

ene = cr_vdwewald(cr,eps);
printf("Edisp (cdcl2 high-p, no damp) : %.10e\n",ene);
