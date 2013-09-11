#! /etc/alternatives/octave -q

addpath("../../../src");

## read from the basin file. Also in critic2 output
x2c = [
       0.92861142170000E+01  -0.46430571085000E+01   0.62554098029316E-15   
       0.00000000000000E+00   0.80420108143658E+01   0.10834687600842E-14
       0.00000000000000E+00   0.00000000000000E+00   0.10215859474400E+02
       ] * .52917720859;
c2x = [
       0.10768766963574E+00   0.62173505052595E-01  -0.13187935992704E-16   
       0.00000000000000E+00   0.12434701010519E+00  -0.13187935992704E-16   
       0.00000000000000E+00   0.00000000000000E+00   0.97887016017194E-01   
       ] / .52917720859;
## from the critic2 output: crystal symmetry operations
op2 = [
       0.000000   -1.000000    0.000000   -0.000000
       1.000000   -1.000000    0.000000    0.000000
       0.000000    0.000000    1.000000    0.333333
       ];
op3 = [
       -1.000000    1.000000    0.000000    0.000000
       -1.000000    0.000000    0.000000    0.000000
       0.000000    0.000000    1.000000   -0.333333
       ];       

## build the underlying crystal structure
cr = cr_read_vasp("POSCAR","POTCAR");
mol = cr_crystalbox(cr,[0.10 0.10 0.0],[1.90 1.90 1.90],2);
rep = representation();
rep = cr_unitcell(cr,rep,[0 0 0],[0 0 1],:,[128 128 128]);
rep = mol_ball(mol,rep,"Si",:,0.4,[100 100 100]);
rep = mol_ball(mol,rep,"O",:,0.25,[255 13 13]);
rep = mol_polyhedron(mol,mol,rep,"Si",{"O"},:,[0 0 128 0 225],[0 0 128],:,:,0.01);

## read and transform the basins.
## the transformations read from critic2's complete list
rep1 = rep_read_basin("quartz-1.basin",:,frgb="",ergb=[255 13 13],:,:,0.0025);
rot = c2x' * op3(:,1:3)' * x2c';
tr = (op3(:,4)' + [1 1 1]) * x2c';
rep1 = rep_transform(rep1,rot,tr);

rep2 = rep_read_basin("quartz-4.basin",:,frgb="",ergb=[255 13 13],:,:,0.0025);
rot = c2x' * op2(:,1:3)' * x2c';
tr = (op2(:,4)' + [1 1 0]) * x2c';
rep2 = rep_transform(rep2,rot,tr);

rep3 = rep_read_basin("quartz-7.basin",:,frgb=[128 128 128 0 0],ergb=[0 0 0],:,:,0.0025);
rot = c2x' * op3(:,1:3)' * x2c';
tr = (op3(:,4)' + [1 1 1]) * x2c';
rep3 = rep_transform(rep3,rot,tr);

## merge all of it
rep = rep_merge(rep,rep1,rep2,rep3);

## write the obj -> surprisingly easy visualization with g3dviewer
## no memory problems at all
rep_write_obj(rep,"quartz_basins.obj");

## orientation and povray
r = [
     0.7660909891 0.1020377576 -0.6345803738 
     -0.6425638199 0.1441552043 -0.7525482178
     0.0146889091 0.9842790365 0.1760020852 
     0.0000000000 0.0000000000 -30.0000000000
     ];
rep = rep_setdefaultscene(rep,r);
save bleh rep

rep_write_pov(rep,"quartz_basins.pov");
system("povray -D -UV +Iquartz_basins.pov +Oquartz_basins.png +W3000 +H3000 +A");

