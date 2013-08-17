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

function cr = cr_newcell(cr0, v1, v2, v3);
% function cr = cr_newcell(cr0, v1, v2, v3)
%
% cr_newcell - change the unit cell.
%
% Required input variables:
% cr0: input crystal.
% v1, v2, v3: crystallographic coordinates of the new cell vectors (1x3).
%
% Output variables:
% cr: output crystal structure.
%

  ## crystal to cartesian
  if (isfield(cr0,"r"))
    r = cr0.r;
    g = r * r';
  else
    if (isfield(cr0,"g"))
      g = cr0.g;
    else
      cc = cos(cr0.b);
      g = cr0.a' * cr0.a;
      g(1,2) = g(2,1) = g(1,2) * cc(3);
      g(1,3) = g(3,1) = g(1,3) * cc(2);
      g(2,3) = g(3,2) = g(2,3) * cc(1);
    endif
    r = chol(g)';
  endif

  mat = [v1; v2; v3];
  minv = inv(mat);
  fvol = det(mat);
  if (abs(fvol) < 1d-6) 
    error("Unit cell with zero volume");
  endif
  
  cr = crystal();
  cr.r = mat' * r;
  cr.g = cr.r * cr.r';
  cr.a = sqrt(diag(cr.g))';
  cr.b(1) = acos(cr.g(2,3) / (cr.a(2)*cr.a(3)));
  cr.b(2) = acos(cr.g(1,3) / (cr.a(1)*cr.a(3)));
  cr.b(3) = acos(cr.g(1,2) / (cr.a(1)*cr.a(2)));
  cr.omega = det(cr.r);
  cr.ntyp = cr0.ntyp;
  cr.ztyp = cr0.ztyp;
  cr.attyp = cr0.attyp;
  if (isfield(cr,"zvaltyp"))
    cr.zvaltyp = cr0.zvaltyp;
  endif

  ## check that fvol is integer
  cr.nat = 0;
  for kk = 0:floor(fvol+1)
    for i = kk:-1:-kk
      for j = kk:-1:-kk
        for k = kk:-1:-kk
          if (abs(i) == kk || abs(j) == kk || abs(k) == kk)
            for l = 1:cr0.nat
              x = cr0.x(l,:) + [i j k];
              xnew = x * minv';
              if (all(xnew > -0.1) && all(xnew < 1.1))
                xnew = xnew - floor(xnew);

                ifound = 0;
                for i0 = 1:cr.nat
                  if (norm(xnew - cr.x(i0,:)) < 1e-10)
                    ifound = 1;
                    break
                  endif
                endfor

                if (ifound == 0) 
                  cr.nat += 1;
                  cr.x(cr.nat,:) = xnew;
                  cr.typ(cr.nat) = cr0.typ(l);
                endif
              endif
            endfor
          endif
        endfor
      endfor
    endfor
  endfor

endfunction
