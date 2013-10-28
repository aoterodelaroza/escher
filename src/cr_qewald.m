% Copyright (c) 2012 Alberto Otero-de-la-Roza and Victor Lua~na
% Adapted from an AOR (2011) routine.
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

function ene = cr_qewald(cr,eps,eta=0)
% function ene = cr_qewald(cr,eps,eta=0)
%
% Electrostatic sum by Ewald method. Outputs energy (Hy) per input cell.
% The value of the point charges is read from the cr.q(1,1:cr.ntyp) array.
% If this is not present, then they are taken from cr.zvaltyp, then from
% cr.ztyp, then an error is raised.
%
% Required input variables:
% {cr}: the crystal description. The following fields are used:
%       a, b, x, typ, qtyp (or zvaltyp or zval).
% {eps}: upper bound for error in real/reciprocal energy, Hy
%
% Optional input variables (all have default values):
% {eta}: if given, non-standard Ewald convergence coefficient
%
% Output:
% {ene}: electrostatic energy per unit cell (Hy).
%

  a = cr.a;
  b = cr.b;
  x = cr.x;
  typ = cr.typ;
  if (isfield(cr,"qtyp") && !isempty(cr.typ))
    q = cr.qtyp;
  elseif (isfield(cr,"zvaltyp") && !isempty(cr.zvaltyp))
    q = cr.zvaltyp;
  elseif (isfield(cr,"ztyp") && !isempty(cr.ztyp))
    q = cr.ztyp;
  else
    error("cr_qewald: charges vector not found.")
  endif

  SGROW = 1.4;
  EPSCUT = 1d-5;

  # find indices for reordering
  irmin = find(a == min(a))(1);
  cb = cos(b);
  G = [ a(1)*a(1),       a(1)*a(2)*cb(3), a(1)*a(3)*cb(2);
        a(1)*a(2)*cb(3), a(2)*a(2),       a(2)*a(3)*cb(1);
        a(1)*a(3)*cb(2), a(2)*a(3)*cb(1), a(3)*a(3) ];
  omega = sqrt(det(G));
  Gr = inv(G);
  ar = sqrt(diag(Gr))';
  ihmin = find(ar == min(ar));
  for i = length(ihmin):-1:1
    if (ihmin(i) != irmin)
      ihmin = ihmin(i);
      break
    endif
  endfor
  for i = 1:3
    if (i != irmin && i != ihmin)
      oidx = i;
      break
    endif
  endfor

  # permute
  if (irmin == 1 && ihmin == 3)
    p = eye(3);
  elseif (irmin == 1 && ihmin == 2)
    p = [1 0 0
         0 0 1
         0 1 0];
  elseif (irmin == 2 && ihmin == 1)
    p = [0 1 0
         1 0 0
         0 0 1] * [1 0 0
                   0 0 1
                   0 1 0];
  elseif (irmin == 2 && ihmin == 3)
    p = [0 1 0
         1 0 0
         0 0 1];
  elseif (irmin == 3 && ihmin == 1)
    p = [0 0 1
         0 1 0
         1 0 0];
  elseif (irmin == 3 && ihmin == 2)
    p = [1 0 0
         0 0 1
         0 1 0] * [0 1 0
                   1 0 0
                   0 0 1];
  endif

  a *= p;
  b *= p;
  x *= p;

  # pre-calc.
  sqpi = sqrt(pi);
  nat = size(x,1);

  qsum = sum(q(typ(1:nat)));
  q2sum = sum(q(typ(1:nat)).^2);

  # calculate cutoffs
  if (eta == 0)
    eta = sqrt(omega / pi / a(2) / sin(b(3)));
  endif
  ##printf("eta = %.15f\n",eta)

  # real space cutoff
  rcut1 = 1;
  rcut2 = 2 / SGROW;
  do
    rcut2 = rcut2 * SGROW;
    err_real = pi * nat^2 * q2sum / omega * eta^2 * erfc(rcut2 / eta);
  until(err_real < eps)
  while(rcut2-rcut1 > EPSCUT)
    rcut = 0.5*(rcut1+rcut2);
    err_real = pi * nat^2 * q2sum / omega * eta^2 * erfc(rcut / eta);
    if (err_real > eps)
      rcut1 = rcut;
    else
      rcut2 = rcut;
    endif
  endwhile
  rcut = 0.5*(rcut1+rcut2);
  ##printf("rcut = %.15f\n",rcut)

  # real space cells to explore
  lrmax = zeros(1,3);
  lrmax(1) = a(2) * a(3) * sin(b(1));
  lrmax(2) = a(1) * a(3) * sin(b(2));
  lrmax(3) = a(1) * a(2) * sin(b(3));
  lrmax = floor(rcut * lrmax / omega) + 1;
  ##printf("lrmax = %d %d %d\n",lrmax)
  
  # reciprocal space cutoff
  hcut1 = 1;
  hcut2 = 2 / SGROW;
  do
    hcut2 = hcut2 * SGROW;
    err_rec = nat^2 * q2sum / sqpi / eta * erfc(eta * hcut2 / 2);
  until(err_rec < eps)
  while(hcut2-hcut1 > EPSCUT)
    hcut = 0.5*(hcut1+hcut2);
    err_rec = nat^2 * q2sum / sqpi / eta * erfc(eta * hcut / 2);
    if (err_rec > eps)
      hcut1 = hcut;
    else
      hcut2 = hcut;
    endif
  endwhile
  hcut = 0.5*(hcut1+hcut2);
  ##printf("hcut = %.15f\n",hcut)

  # reciprocal space cells to explore
  lhmax = floor(a / (2*pi) * hcut) + 1;
  ##printf("lhmax = %d %d %d\n",lhmax)

  # real sum
  rcut2 = rcut * rcut;
  sum_real = 0;
  for i1 = -lrmax(1):lrmax(1)
    for i2 = -lrmax(2):lrmax(2)
      for i3 = -lrmax(3):lrmax(3)
        L = [i1 i2 i3];

        for i = 1:nat
          for j = 1:nat
            px = x(i,:) - x(j,:) - L;
            d2 = (px * G * px');
            if (d2 < 1d-12 || d2 > rcut2)
              continue
            endif
            d = sqrt(d2) / eta;
            sum_real += q(typ(i)) * q(typ(j)) * erfc(d) / d;
          endfor
        endfor

      endfor
    endfor
  endfor

  sum_real /= 2 * eta;
  ##printf("sum_real = %.15f\n",sum_real);

  # reciprocal sum
  sum_rec = 0;
  for i1 = -lhmax(1):lhmax(1)
    for i2 = -lhmax(2):lhmax(2)
      for i3 = -lhmax(3):lhmax(3)
        h = 2 * pi * [i1 i2 i3];
        dh = sqrt(h * Gr * h');
        if (dh < 1d-12 || dh > hcut)
          continue
        endif

        bb = 0.5 * dh * eta;

        sfac_s = 0; sfac_c = 0;
        for i = 1:nat
          sfac_s += q(typ(i)) * sin(h * x(i,:)');
          sfac_c += q(typ(i)) * cos(h * x(i,:)');
        endfor
        sfacp = sfac_s^2 + sfac_c^2;

        sum_rec += sfacp / dh^2 * exp(-bb^2);
      endfor
    endfor
  endfor
  sum_rec *= 2 * pi / omega;
  ##printf("sum_rec = %.15f\n",sum_rec);

  # h = 0 term
  sum0 = - q2sum / sqpi / eta;
  ##printf("sum_0 = %.15f\n",sum0);

  # compensating background charge term
  sum_back = - qsum^2 * eta^2 * pi / omega / 2;

  ene = sum_real + sum_rec + sum0 + sum_back;

endfunction

