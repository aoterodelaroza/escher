mol_dbstart(0);

for n = 3:5;
dPN = 1.58;
dPX = 2.00;
aXPX = 103.0;
beta = (pi-pi/n)/2;
dring = dPN / (2*cos(beta));
for i = 1 : 2*n
   phi = (i-1) * pi / n;
   xyz = [dring*sin(phi); dring*cos(phi); 0];
   if (mod(i,2) == 1)
      if (i==1)
         frg1 = mol_addatom("P",xyz,[],1);
      else
         frg1 = mol_addatom("P",xyz,frg1);
      endif
      r = dring + dPX * cos(aXPX/2 * pi/180);
      z = dPX * sin(aXPX/2 * pi/180);
      frg1 = mol_addatom("Cl",[r*sin(phi);r*cos(phi);+z],frg1);
      frg1 = mol_addatom("Cl",[r*sin(phi);r*cos(phi);-z],frg1);
   else
      frg1 = mol_addatom("N",xyz,frg1);
   endif
endfor

# planar cycle:
file = sprintf("pnx2_%02d_plane.xyz", n);
mol_writexyz(frg1, file, "anxyz");
endfor
