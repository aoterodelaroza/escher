#! /usr/bin/octave -q

addpath("../../../src");
cr = cr_read_espresso("serine.scf.out");
printf("Title : %s\n",cr.name);
printf("Cell lengths (bohr): %f %f %f\n",cr.a);
printf("Cell angles (deg): %f %f %f\n",cr.b * 180 / pi);
printf("Volume (bohr^3): %f\n",cr.omega);
printf("Atomic positions: \n");
for i = 1:cr.nat
  printf("%s %f %f %f\n",cr.attyp{cr.typ(i)},cr.x(i,:));
endfor


