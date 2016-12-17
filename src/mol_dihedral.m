% Copyright (C) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function [dh x1p x2p] = mol_dihedral(mol, at1, at2, at3, at4, extend=[0 0 0])
% function [dh x1p x2p] = mol_dihedral(mol, at1, at2, at3, at4, extend=[0 0 0])
%
% mol_dihedral - returns the dihedral angle (in degrees) between four
% atoms in the molecule. at2-at3 is the central bond. If x1p and x2p
% are present, return encompassing quadrangles for each of the two
% planes forming the dihedral. This is useful in combination with
% rep_polygon. Extend the quadrangles' axes in the x-direction and the
% two y-directions by a distance given by the extend argument
%
% Required input variables:
% mol: structure containing the molecule.
% at1, at2, at3, at4: return the at1-at2-at3-at4 dihedral angle (in degrees).
% extend: extend the dihedral planes by this distance in the x (first element),
% and the two y directions (second and third).
%
% Authors: VLC Victor Lua~na .......... <victor@fluor.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <aoterodelaroza@gmail.com>
% Created: Nov 2012

  b1 = mol.atxyz(:,at2)-mol.atxyz(:,at1); 
  b2 = mol.atxyz(:,at3)-mol.atxyz(:,at2); 
  b3 = mol.atxyz(:,at4)-mol.atxyz(:,at3); 
  b12 = cross(b1, b2);
  b23 = cross(b2, b3);
  dh = atan2(norm(b2) * dot(b1, b23), dot(b12, b23)) * (180/pi);

  if (nargout() > 1)
    x0 = mol.atxyz(:,at2)';
    x1 = (mol.atxyz(:,at3)-mol.atxyz(:,at2))';
    x2 = (mol.atxyz(:,at1)-mol.atxyz(:,at2))';
    x3 = (mol.atxyz(:,at4)-mol.atxyz(:,at2))';

    ux = x1 / norm(x1);
    uy1 = x2 - (x2 * ux') * ux;
    uy1 = uy1 / norm(uy1);
    uy2 = x3 - (x3 * ux') * ux;
    uy2 = uy2 / norm(uy2);
    
    xcoord = [x1 * ux', x2 * ux', x3 * ux'];
    ycoord1 = [x1 * uy1', x2 * uy1', x3 * uy1'];
    ycoord2 = [x1 * uy2', x2 * uy2', x3 * uy2'];
    uxmin = min(xcoord) - extend(1);
    uymin1 = min(ycoord1) - extend(2);
    uymin2 = min(ycoord2) - extend(3);
    lx = range(xcoord) + 2 * extend(1);
    ly1 = range(ycoord1) + 2 * extend(2);
    ly2 = range(ycoord2) + 2 * extend(3);

    x0 = (mol.atxyz(:,at2)' + uxmin * ux + uymin1 * uy1);
    x1p(1,:) = x0;
    x1p(2,:) = (x0 + lx * ux);
    x1p(3,:) = (x0 + lx * ux + ly1 * uy1);
    x1p(4,:) = (x0 + ly1 * uy1);

    x0 = (mol.atxyz(:,at2)' + uxmin * ux + uymin2 * uy2);
    x2p(1,:) = x0;
    x2p(2,:) = (x0 + lx * ux);
    x2p(3,:) = (x0 + lx * ux + ly2 * uy2);
    x2p(4,:) = (x0 + ly2 * uy2);
  endif

endfunction
