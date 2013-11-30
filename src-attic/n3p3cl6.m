mol_dbstart(0);

n = 3;
alfa = pi/n;
beta = (pi-alfa)/2;
dpn = 1.1;
dpcl = 1.2;
frg = mol_addatom("S",[0.0, 0.0, 0.0]',[],1);
frg = mol_addatom("C",[0.0, 0.0, 1.1]',frg);
xyz = [0.5, 0.0, 1.6]';
c20 = op_rotz(20);
xyz = c20 * xyz;
frg = mol_addatom("H",xyz,frg);
c180 = op_rotz(180);
xyz = c180 * xyz;
frg = mol_addatom("H",xyz,frg);
cx = op_rotx(45);
frg = mol_transform(op_rotx(45), frg);

c6 = op_rotz(60);
xyz = [1.2, 0.0, -2.3]';
for i = 1:6
   frg = mol_addatom("C",xyz,frg);
   xyz = c6 * xyz;
endfor
xyz = [2.2, 0.0, -2.3]';
for i = 1:6
   frg = mol_addatom("H",xyz,frg);
   xyz = c6 * xyz;
endfor
file = sprintf("planar_%02d.xyz", n);
mol_writexyz(frg);
