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

function cr = cr_insertvacuum(cr0, ix, r0, chop="", center="");
% function cr = cr_insertvacuum(cr0, ix, r0, chop="", center="")
%
% cr_insertvacuum - insert vacuum in one of the crystallographic directions. 
% This is useful for surface design.
%
% Required input variables:
% cr0: input crystal.
% ix: vacuum direction (1=x, 2=y, 3=z).
% r0: amount of vacuum, in bohr.
%
% Optional input variables:
% chop: if a range [x0 x1] is given, use only the fraction of the
%       cr atoms which are between x0 and x1 along the ix axis before
%       inserting the vacuum.
%       This is useful when building centrosymmetric surfaces.
% center: the origin of the new cell is displaced half the vacuum length (r0) 
%         in the ix direction. If the cell of cr0 had an inversion center at
%         the origin, then the new cell does as well. Recommended if you 
%         are building a centrosymmetric surface, since solid-state
%         codes like this.
%
% Output variables:
% cr: output crystal structure.
%

  if (!any(ix == [1 2 3]))
    error("ix must be one of 1 (x), 2 (y), or 3 (z)");
  endif

  ## crystal to cartesian
  r = cr0.r;
  g = cr0.g;

  ## insert the vacuum 
  aa0 = norm(r(ix,:));
  aa = aa0 + r0;
  r(ix,:) = r(ix,:) / aa0 * aa;

  ## create the new crystal and renormalize all atomic coordinates
  cr = cr0;
  cr.r = r;
  cr.g = cr.r * cr.r';
  cr.a = sqrt(diag(cr.g))';
  cr.b(1) = acos(cr.g(2,3) / (cr.a(2)*cr.a(3)));
  cr.b(2) = acos(cr.g(1,3) / (cr.a(1)*cr.a(3)));
  cr.b(3) = acos(cr.g(1,2) / (cr.a(1)*cr.a(2)));
  cr.omega = det(cr.r);
  cr.nat = 0;
  cr.x = [];
  cr.typ = [];
  for i = 1:cr0.nat
    if (!isempty(chop) && (cr0.x(i,ix) < chop(1) || cr0.x(i,ix) > chop(2))) 
      continue
    endif
    cr.nat += 1;
    cr.x(cr.nat,:) = cr0.x(i,:);
    cr.x(cr.nat,ix) = cr0.x(i,ix) * aa0 / aa;
    cr.typ(cr.nat) = cr0.typ(i);
  endfor
  if (!isempty(center))
    tr = [0 0 0];
    tr(ix) = -r0/aa/2;
    cr = cr_moveorigin(cr,tr);
  endif

endfunction
