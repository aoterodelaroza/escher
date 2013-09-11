#! /etc/alternatives/octave -q

addpath("../../../src");
cr = cr_read_vasp("POSCAR","POTCAR",1);

mol = cr_crystalbox(cr,x0=[-1.05 -1.05 -1.05],x1=[2.05 2.05 2.05]);
mol_writexyz(mol,"agtao3.xyz");

