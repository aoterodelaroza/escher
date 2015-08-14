mol_dbstart(0);

frg1 = mol_addatom("S",[0.0, 0.0, 0.0]',[],1);
frg1 = mol_addatom("C",[0.0, 0.0, 1.1]',frg1);
xyz = [0.5, 0.0, 1.6]';
c20 = op_rotz(20);
xyz = c20 * xyz;
frg1 = mol_addatom("H",xyz,frg1);
c180 = op_rotz(180);
xyz = c180 * xyz;
frg1 = mol_addatom("H",xyz,frg1);
cx = op_rotx(45);
frg1 = mol_transform(op_rotx(45), frg1);

c6 = op_rotz(60);
xyz = [1.2, 0.0, -2.3]';
for i = 1:6
   frg1 = mol_addatom("C",xyz,frg1);
   xyz = c6 * xyz;
endfor
xyz = [2.2, 0.0, -2.3]';
for i = 1:6
   frg1 = mol_addatom("H",xyz,frg1);
   xyz = c6 * xyz;
endfor
mol_writexyz(frg1);

printf("C(2)-H(3): %.4f\n", norm(frg1.atxyz(:,2)-frg1.atxyz(:,3)));

