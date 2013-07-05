% Copyright (C) 2011--12 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [mol] = mol_adduct(mol1, mol2, iat1=-1, iat2=-1, rtp1=[0,0,0], iang1=[0,0,0], order1='xyz', rtp2=[0,0,0], iang2=[0,0,0], order2='xyz', LOG=1)
% function [mol] = mol_adduct(mol1, mol2, iat1=-1, iat2=-1, rtp1=[0,0,0], iang1=[0,0,0], order1='xyz', rtp2=[0,0,0], iang2=[0,0,0], order2='xyz', LOG=1)
%
% mol_adduct - create the adduct of mol1 and mol2.
%
% Required input variables:
% {mol1,mol2}: molecules to join.
%
% Optional input variables (all have default values):
% {iat1,iat2}: reference points for mol1 and mol2. It can be:
%     -1) no ref. point = don't move the molecule (before applying the
%         next transformations.
%      0) the center of mass of the input molecule.
%     >0) number of the atom that will be moved to the origin.
% {rtp1,rtp2}: (r,theta,phi) spherical coordinates of the final position
%     to be occupied by the reference point in mol1 and mol2.
% {iang1,iang2}: (alpha,beta,gamma) Euler-like angles that describe the
%     orientation of mol1 and mol2, respectively, in the final adduct.
% {order1,order2}: order to apply the orientational rotation of the mol1
%     and mol2 fragments. An order of "xyz" would correspond to multiply
%     the column coordinates by:
%     vector = op_rotz(gamma) * op_roty(beta) * op_rotx(alpha) * vector.
% {LOG}: print the final result if LOG>0.
%
% Required output variables:
% {mol}: molecular database for the adduct.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

n1 = length(mol1.atmass);
if (iat1 == 0)
   # Center of mass of mol1:
   cm1 = mol1.atxyz * mol1.atmass' / sum(mol1.atmass);
   # Move mol1 to place cm1 at the origin:
   x1 = mol1.atxyz - cm1 * ones(1,n1);
elseif (iat1 > 0)
   # Move mol1 to place iat1 at the origin:
   if (iat1 > n1)
      error('mol_adduct: iat1 atom is beyond limits!')
   endif
   x1 = mol1.atxyz - mol1.atxyz(3,iat1) * ones(1,n1);
elseif (iat1 < 0)
   # Don't move anything
   x1 = mol1.atxyz;
endif

n2 = length(mol2.atmass);
if (iat2 == 0)
   # Center of mass of mol2:
   cm2 = mol2.atxyz * mol2.atmass' / sum(mol2.atmass);
   # Move mol2 to place cm2 at the origin:
   x2 = mol2.atxyz - cm2 * ones(1,n2);
elseif (iat2 > 0)
   # Move mol2 to place iat2 at the origin:
   if (iat2 > n2)
      error('mol_adduct: iat2 atom is beyond limits!')
   endif
   x2 = mol2.atxyz - mol2.atxyz(3,iat2) * ones(1,n2);
elseif (iat2 < 0)
   # Don't move anything
   x2 = mol2.atxyz;
endif

# Apply internal rotations to mol1:
if (isequal(order1,'xyz'))
   op = op_rotz(iang1(3)) * op_roty(iang1(2)) * op_rotx(iang1(1));
elseif (isequal(order1,'xzy'))
   op = op_roty(iang1(3)) * op_rotz(iang1(2)) * op_rotx(iang1(1));
elseif (isequal(order1,'yxz'))
   op = op_rotz(iang1(3)) * op_rotx(iang1(2)) * op_roty(iang1(1));
elseif (isequal(order1,'yzx'))
   op = op_rotx(iang1(3)) * op_rotz(iang1(2)) * op_roty(iang1(1));
elseif (isequal(order1,'zxy'))
   op = op_roty(iang1(3)) * op_rotx(iang1(2)) * op_rotz(iang1(1));
elseif (isequal(order1,'zyx'))
   op = op_rotx(iang1(3)) * op_roty(iang1(2)) * op_rotz(iang1(1));
else
   error('mol_adduct: internal rotation order1 not defined!');
endif
x1 = op(1:3,1:3)*x1;

# Move reference point in mol1 to its final position:
r = rtp1(1); theta = rtp1(2); phi = rtp1(3);
x1t = [r*sind(theta)*cosd(phi), r*sind(theta)*sind(phi), r*cosd(theta)]';
x1 = x1 + x1t * ones(1,n1);

# Apply internal rotations to mol2:
if (isequal(order2,'xyz'))
   op = op_rotz(iang2(3)) * op_roty(iang2(2)) * op_rotx(iang2(1));
elseif (isequal(order2,'xzy'))
   op = op_roty(iang2(3)) * op_rotz(iang2(2)) * op_rotx(iang2(1));
elseif (isequal(order2,'yxz'))
   op = op_rotz(iang2(3)) * op_rotx(iang2(2)) * op_roty(iang2(1));
elseif (isequal(order2,'yzx'))
   op = op_rotx(iang2(3)) * op_rotz(iang2(2)) * op_roty(iang2(1));
elseif (isequal(order2,'zxy'))
   op = op_roty(iang2(3)) * op_rotx(iang2(2)) * op_rotz(iang2(1));
elseif (isequal(order2,'zyx'))
   op = op_rotx(iang2(3)) * op_roty(iang2(2)) * op_rotz(iang2(1));
else
   error('mol_adduct: internal rotation order2 not defined!');
endif
x2 = op(1:3,1:3)*x2;

# Move reference point in mol2 to its final position:
r = rtp2(1); theta = rtp2(2); phi = rtp2(3);
x2t = [r*sind(theta)*cosd(phi), r*sind(theta)*sind(phi), r*cosd(theta)]';
x2 = x2 + x2t * ones(1,n2);

# Form, finally, the adduct:
mol.name = strcat('Adduct: ', mol1.name, ' + ', mol2.name);
new = 0;
for i = 1 : n1
   new++;
   mol.atname{new} = mol1.atname{i};
   mol.atnumber(new) = mol1.atnumber(i);
   mol.atmass(new) = mol1.atmass(i);
   mol.atxyz(1:3,new) = x1(1:3,i);
endfor
for i = 1 : n2
   new++;
   mol.atname{new} = mol2.atname{i};
   mol.atnumber(new) = mol2.atnumber(i);
   mol.atmass(new) = mol2.atmass(i);
   mol.atxyz(1:3,new) = x2(1:3,i);
endfor

if (LOG > 0)
   printf("Adduct from %s and %s\n", mol1.name, mol2.name);
   printf("Reference for mol%d: %d\n", 1, iat1);
   printf("Orientation (%s) for mol%d: %.3f %.3f %.3f\n", order1, 1, iang1);
   printf("Position for mol%d: %.3f %.3f %.3f\n", 1, rtp1);
   printf("Reference for mol%d: %d\n", 2, iat2);
   printf("Orientation (%s) for mol%d: %.3f %.3f %.3f\n", order2, 2, iang2);
   printf("Position for mol%d: %.3f %.3f %.3f\n", 2, rtp2);
endif

endfunction
