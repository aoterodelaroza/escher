#! /etc/alternatives/octave -q

addpath("../../../src");

# read and extract a 11-ang-radius sphere around the
# position of the methane.
cr = cr_read_espresso("clathrate.scf.out",0);
mol = cr_spherebox(cr,[0 0 0.25],11);

# keep the hydrogens directly attached to something
trim = [];
for i = 1:length(mol.atnumber)
  if (mol.atnumber(i) == 1)
    for j = 1:length(mol.atnumber)
      if (mol.atnumber(j) > 1 && (norm(mol.atxyz(:,i)-mol.atxyz(:,j)) < 1.2))
        trim = [trim i];
        break
      endif
    endfor
  else
    trim = [trim i];
  endif
endfor
mol = mol_getfragment(mol,trim);


# create a reprsentation where the clathrate cage is evident,
# with covalent bonds in blue and hydrogen bonds in green.
rep = mol_ball(mol,:,:,:,-0.3);
rep = mol_stick(mol,rep,"H","O",[0.9 1.1],:,:,[0 255 0]);
rep = mol_stick(mol,rep,"H","O",[1.3 1.8],:,:,[0 0 255]);
rep = mol_stick(mol,rep,"H","C",[0.8 1.2],:,:,[255 0 0]);
rep = rep_setdefaultscene(rep);
rep = rep_setbgcolor(rep,[0 0 0]);

# geomview output
rep_write_off(rep,"clathrate.off");

