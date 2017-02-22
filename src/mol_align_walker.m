% Copyright (C) 2011 Victor Lua~na and Alberto Otero-de-la-Roza
%
% This octave routine is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or (at
% your option) any later version. See <http://www.gnu.org/licenses/>.
%
% The routine distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for
% more details.

function [molout u rmsd] = mol_align_walker(mol0, mol1)
% function [molout u rmsd] = mol_align_walker(mol0, mol1)
%
% mol_align_walker - rotate the mol0 molecule in order to minimize
% the rmsd of the atomic positions compared to the mol1 molecule. The
% two molecules must have the same number of atoms and the atoms must
% be in the same order. The output molecule is the mol0 molecule,
% rotated and translated to the mol1 center of mass.  The rotation
% matrix u, and the RMSD of the atomic positions (rmsd) are also
% given.  This routine uses Walker's algorithm based on quaternion
% algebra. (Walker et al., CVGIP-Imag. Understan. 54 (1991) 358.)
%
% Input variables:
% mol0: rotated molecule
% mol1: target molecule.
% allow_inversion: if 1, permit the application of an inversion
% operation to match the two molecules. Warning: if mol0 is chiral, it
% is converted to its enantiomer.
%
% Output variables:
% molout: rotated and translated mol0.
% u: rotation matrix.
% rmsd: 
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: June 2011

  if (mol0.nat != mol1.nat)
    error("inconsistent number of atoms in mol0 and mol1")
  endif
  xcm0 = mol_cmass(mol0);
  xcm1 = mol_cmass(mol1);
  mol0.atxyz = mol0.atxyz - xcm0 * ones(1,mol0.nat);
  mol1.atxyz = mol1.atxyz - xcm1 * ones(1,mol1.nat);

  ## the W and Q functions (eqs. 15 and 16)
  wmat = @(x) [
               x(4), x(3), -x(2), x(1),
               -x(3), x(4), x(1), x(2),
               x(2), -x(1), x(4), x(3),
               -x(1), -x(2), -x(3), x(4),
          ];
  qmat = @(x) [
               x(4), -x(3), x(2), x(1),
               x(3), x(4), -x(1), x(2),
               -x(2), x(1), x(4), x(3),
               -x(1), -x(2), -x(3), x(4),
          ];

  
  ## The c1, c2, and c3 quaternions (eqs. 35-37)
  n = mol0.nat;
  c2 = n / 2;
  c1 = c3 = zeros(4);
  for i = 1:n
    x = [mol0.atxyz(:,i)' 0];
    w = wmat(x);
    x = [mol1.atxyz(:,i)' 0];
    q = qmat(x);
    c1 = c1 - (q' * w);
    c3 = c3 + (w - q);
  endfor

  ## the a matrix (eq. 47), diagonalization, and rotation matrix
  a = (c3' * c3) .* c2 - c1;
  [v, p] = eig(a);
  x = v(:,find(diag(p) == max(diag(p))));
  qf = wmat(x)' * qmat(x);
  
  ## output quantities
  u = qf(1:3,1:3);
  atxyz0 = u * mol0.atxyz;
  rmsd = sqrt(sum(sum((atxyz0 - mol1.atxyz).^2)) / mol0.nat);
  molout = mol0;
  molout.atxyz = atxyz0 + xcm1 * ones(1,n);

endfunction
