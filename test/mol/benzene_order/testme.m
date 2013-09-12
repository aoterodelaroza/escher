#! /usr/bin/octave -q

addpath("../../../src");

# This example shows how to use mol_burst to detect the molecular
# fragments in a mol and how to reorder the atoms. The problem is: given
# benzene clusters up to n molecules, with n=3 to 6, separate the xyz
# into molecules (mol_burst), reorder in descending atomic number
# (mol_order) and then repack them again in a cluster (mol_merge). 

for n = 3:6
  smol = mol_burst(mol_readxyz(sprintf("b0%d.xyz",n)));
  mol = molecule();
  for i = 1:length(smol)
    mol = mol_merge(mol,mol_order(smol{i},'atnumber',@fliplr));
  endfor
  mol_writexyz(mol,sprintf("b0%d_ordered.xyz",n));
endfor





