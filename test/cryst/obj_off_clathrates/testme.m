#! /usr/bin/octave -q

addpath("../../../src");

## load the crystal information and get the scale
load("crystal_info.dat","cr_all");
emin = 0;
emax = -1e30;
for i = 1:7
  for j = 1:7
    emin = min(cr_all{i,j}.energy,emin);
    emax = max(cr_all{i,j}.energy,emax);
  endfor
endfor
r2k = 627.52 / 2;
scale = ([emin emax] - emin) * r2k;

## Extract the cage. Discard all hydrogens not directly bonded to an oxygen. 
cr = cr_all{1,1};
mol = cr_spherebox(cr,[0 0 0.25],[3 11]);
trim = [];
for i = 1:length(mol.atnumber)
  if (mol.atnumber(i) == 1)
    for j = 1:length(mol.atnumber)
      if (mol.atnumber(j) > 1 && norm(mol.atxyz(:,i)-mol.atxyz(:,j)) < 1.2)
        trim = [trim i];
        break
      endif
    endfor
  else
    trim = [trim i];
  endif
endfor
cage = mol_getfragment(mol,trim);

## Create the representation for the cage and set up the scene
rcageh = mol_ball(cage,:,:,:,-0.3);
rcageh = mol_stick(cage,rcageh,"H","O",[0.9 1.1],:,:,[0 255 0]);
rcageh = mol_stick(cage,rcageh,"H","O",[1.3 1.8],:,:,[0 0 255]);
rcageh = rep_setdefaultscene(rcageh,0);
rcageh = rep_setbgcolor(rcageh,[0 0 0]);

## accumulate all CO2 molecules in rmol
rmol = representation();
plist = [0 30 60 90 120 150 180];
tlist = [0 30 60 90 120 150 180];
for ip = 1:length(plist)
  for it = 1:length(tlist)
    ## energy to color
    e = (cr_all{ip,it}.energy - emin)*r2k;
    se = (e - scale(1))/scale(2);
    rgb = max(round([sqrt(se), se**3, sin(2*pi*se)] * 255),0);

    ## position, add the C atom only once.
    cr = cr_all{ip,it};
    mol = cr_spherebox(cr,[0 0 0.25],3);
    if (ip == 1 && it == 1)
      rmol = mol_ball(mol,rmol,"C",0,-0.3);
    endif
    rmol = mol_ball(mol,rmol,"O",0,-0.3,rgb);
  endfor
endfor

## join the cage and the CO2s and plot
rep = rep_merge(rcageh,rmol);
rep_write_off(rep,"clathrate_co2_h_rot.off");
rep_write_obj(rep,"clathrate_co2_h_rot.obj");
