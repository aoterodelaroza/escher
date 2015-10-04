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

function ene = cr_vdwewald(cr,eps,idamp=0,iparam=[])
% function ene = cr_vdwewald(cr,eps,idamp=0,iparam=[])
%
% Dispersion sum (R^-6) by Ewald method. Outputs energy (Hy) per input cell. 
% The interaction coefficients (c6) are read from cr.c6typ(cr.ntyp,cr.ntyp).
% Several short-range damping schemes can be used:
% * idamp=0, no damping.
% * idamp=1, fermi-type f(d) = 1 / (1+exp(-iparam(1)*(d/iparam(2)/(r1+r2)-1)))
% * idamp=2, becke-johnson f(d) = (d^6 / (d^6 + (iparam(1)*r1+iparam(2)+iparam(1)*r2+iparam(2))^6));
% where r1 and r2 are the van der waals radii for atoms 1 and 2 (that
% is, cr.rvdwtyp(typ(:))), and d is the interatomic distance. Note that the interpretation
% of iparam depends on the type of damping.
%
% If anything other than idamp=0 is used, then the field cr.rvdwtyp(cr.ntyp,cr.ntyp),
% which contains the van der Waals radii is used. The iparam array contains
% the parameters of the damping, and its interpretation depends on the value of idamp.
%
% Required input variables:
% {cr}: the crystal description. The following fields are used:
%       a, b, x, typ, c6typ, cr.rvdwtyp (optional).
% {eps}: upper bound for error in real/reciprocal energy, Hy
%
% Optional input variables (all have default values):
% {idamp}: type of damping function (see above).
% {iparam}: array of damping parameters (see above).
%
% Output:
% {ene}: dispersion energy per unit cell (Hy).
%

  a = cr.a;
  b = cr.b;
  x = cr.x;
  typ = cr.typ;
  c6 = cr.c6typ;
  if (isfield(cr,"rvdwtyp")) 
    rvdw = cr.rvdwtyp;
  else
    rvdw = zeros(cr.ntyp);
  endif

  SGROW = 1.4;
  EPSCUT = 1d-5;

  if (idamp == 1)
    ## fermi type damping (tkatchenko,...)
    fdamp = @(d,r1,r2)(1 ./ (1+exp(-iparam(1)*(d/iparam(2)/(r1+r2)-1))));
  elseif (idamp == 2)
    ## becke-johnson damping
    fdamp = @(d,r1,r2)(d.^6 ./ (d.^6 + (iparam(1)*r1+iparam(2)+iparam(1)*r2+iparam(2)).^6));
  endif

  # find indices for reordering
  irmin = find(a == min(a))(1);
  cb = cos(b);
  G = [ a(1)*a(1),       a(1)*a(2)*cb(3), a(1)*a(3)*cb(2);
        a(1)*a(2)*cb(3), a(2)*a(2),       a(2)*a(3)*cb(1);
        a(1)*a(3)*cb(2), a(2)*a(3)*cb(1), a(3)*a(3) ];
  omega = sqrt(det(G));
  Gr = inv(G);
  ar = sqrt(diag(Gr))';
  ihmin = find(ar < min(ar)+1d-12);
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
  else
    p = eye(3);
  endif

  a *= p;
  b *= p;
  x *= p;

  # pre-calc.
  sqpi = sqrt(pi);
  nat = size(x,1);
  mtyp = max(typ);

  c6sum = sum(sum(c6(typ(1:nat),typ(1:nat))));

  # calculate cutoffs
  eta = sqrt(omega / pi / a(2) / sin(b(3)));
  ##printf("eta = %.15f\n",eta)

  # real space cutoff
  rcut1 = 1;
  rcut2 = 2 / SGROW;
  do
    rcut2 = rcut2 * SGROW;
    err_real = pi^(3/2) * c6sum * eta / omega * ...
        (1/rcut2^4 + 1/rcut2^2/eta^2 + 0.5/eta^4) * ...
        erfc(rcut2/eta);
  until(err_real < eps)
  while(rcut2-rcut1 > EPSCUT)
    rcut = 0.5*(rcut1+rcut2);
    err_real = pi^(3/2) * c6sum * eta / omega * ...
        (1/rcut^4 + 1/rcut^2/eta^2 + 0.5/eta^4) * ...
        erfc(rcut/eta);
    if (err_real > eps)
      rcut1 = rcut;
    else
      rcut2 = rcut;
    endif
  endwhile
  rcut = 0.5*(rcut1+rcut2);
  if (idamp > 0)
    ## real space cutoff
    rcut = 0;
    for i = 1:length(rvdw)
      for j = i:length(rvdw)
        rcut = max(3*(rvdw(i)+rvdw(j)),rcut);
      endfor
    endfor
  endif
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
    err_rec = c6sum / 6 / sqpi / eta^6 * ...
        (hcut2 * eta * exp(-(hcut2 * eta / 2)^2) + ...
         sqpi * erfc(hcut2 * eta / 2));
  until(err_rec < eps)
  while(hcut2-hcut1 > EPSCUT)
    hcut = 0.5*(hcut1+hcut2);
    err_rec = c6sum / 6 / sqpi / eta^6 * ...
        (hcut * eta * exp(-(hcut * eta / 2)^2) + ...
         sqpi * erfc(hcut * eta / 2));
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
  sum_damp = 0;
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
            if (idamp > 0)
              sum_damp += c6(typ(i),typ(j)) / d2^3 * ...
                  (1-fdamp(sqrt(d2),rvdw(typ(i)),rvdw(typ(j))));
            endif

            d2 /= eta^2;
            d_2 = 1 / d2;
            sum_real += c6(typ(i),typ(j)) * d_2 * exp(-d2) * (d_2^2 + d_2 + 0.5);
          endfor
        endfor

      endfor
    endfor
  endfor

  sum_real /= 2 * eta^6;
  sum_damp /= 2;
  ##printf("sum_real = %.15f\n",sum_real);
  ##printf("sum_damp = %.15f\n",sum_damp);

  # reciprocal sum
  sum_rec = 0;
  for i1 = -lhmax(1):lhmax(1)
    for i2 = -lhmax(2):lhmax(2)
      for i3 = -lhmax(3):lhmax(3)
        h = 2 * pi * [i1 i2 i3];
        dh = sqrt(h * Gr * h');
        if (dh < 1d-12 || dh > hcut2)
          continue
        endif
        bb = 0.5 * dh * eta;

        sfac = 0;
        for i = 1:nat
          for j = 1:nat
            sfac += c6(typ(i),typ(j)) * cos(h * (x(i,:)-x(j,:))');
          endfor
        endfor
        sum_rec += sfac * dh^3 * (sqpi * erfc(bb) + (0.5 / bb^3 - 1/bb) * exp(-bb^2));

      endfor
    endfor
  endfor
  sum_rec *= pi^(3/2) / 24 / omega;
  ##printf("sum_rec = %.15f\n",sum_rec);

  # h = 0 term
  sum0 = c6sum * pi^(3/2) / omega / 6 / eta^3 - sum(diag(c6)(typ(1:nat))) / 12 / eta^6;
  ##printf("sum_0 = %.15f\n",sum0);

  ene = -(sum_real + sum_rec + sum0 - sum_damp);

endfunction
