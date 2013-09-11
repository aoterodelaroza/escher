#! /etc/alternatives/octave -q

addpath("../../../src");

## NaCl
cr = struct();
cr.a = [5.64056 5.64056 5.64056] / .52917721;
cr.b = [90 90 90] * pi / 180;
cr.nat = 8;
cr.ntyp = 2;
cr.x = [0 0 0
        1/2 1/2 0
        1/2 0 1/2
        0 1/2 1/2
        1/2 1/2 1/2
        1/2 0 0
        0 1/2 0
        0 0 1/2];
cr.typ = [1 1 1 1 2 2 2 2];
cr.ztyp = [11 17];
cr.qtyp = [1 -1];
lambda = 1.5406 ;
npts = 10001;

[t ih pinfo] = cr_xrd(cr,"nacl");

## diamond 
cr = struct();
cr.a = [3.56658 3.56658 3.56658] / .52917721;
cr.b = [90 90 90] * pi / 180;
cr.nat = 8;
cr.ntyp = 1;
cr.x = [1/8 1/8 1/8
        5/8 5/8 1/8
        5/8 1/8 5/8
        1/8 5/8 5/8
        7/8 7/8 7/8
        7/8 3/8 3/8
        3/8 7/8 3/8
        3/8 3/8 7/8];
cr.typ = [1 1 1 1 1 1 1 1];
cr.ztyp = [6];
cr.qtyp = [0];
lambda = 1.5406;
  
[t ih pinfo] = cr_xrd(cr,"diamond");
