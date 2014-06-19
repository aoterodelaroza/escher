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

function [spow S dmax dav delta] = cr_compare(cr1, cr2)
% function [spow S dmax dav delta] = cr_compare(cr1, cr2)
%
% cr_compare - calculate the similarity of two crystal structures using
% popular measures. For S, dmax, dav, and delta, the two crystals can be
% described by any unit cell and the atoms must come in the same order. 
% For spow, this is not necessary. S, dmax, dav, and delta have
% been adapted from the Bilbao crystallographic server.
% 
% Required input variables:
% cr1: first crystal structure.
% cr2: second crystal structure.
%
% Output variables:
% spow - similarity using integrated powder diffraction patterns
% Hofmann and Kuleshova, J. Appl. Cryst. 38 (2005) 861.
% spow = 1/(tmax-tmin) * \int_tmin^tmax 
%      |1/N1 \int_tmin^t0 I1(t) dt - 1/N2 \int_tmin^t0 I2(t) dt| dt0
% Ni = \int_tmin^tmax Ii(t) dt
%
% S - degree of lattice distortion (from the Bilbao crystallographic server)
%     C. Capillas et al. J. Phys.: Condens. Matter 19 (2007) 275203.
% S = 1/3 * sqrt(eta_1 + eta_2 + eta_3)
% eta_i = eigenvalues of eta (finite Lagrangian strain tensor)
% eta = 1/2 (e + e^t + e^t * e)
% e = R2 * R1^-1 - I
% R = standard root tensor (cr.r)
%
% dmax - maximum distance (in bohr)
% dmax = max(|x1 - x2|)
%
% dav = average distance (in bohr)
% dav = mean(|x1 - x2|)
%
% delta - measure of similarity
% Bergerhoff et al. Acta Cryst. (1999) B55 147-156.
% delta = (sqrt(2) * deltac + 1) * deltad - 1
% deltac = sum(m_i * |x1_i - x2_i|^2)/sum(m_i)
% deltad = (b1/a1)*(c1/a1)/(b2/a2)/(c2/a2)
%
  
  if (cr1.nat == cr2.nat) 
    ## calculate S
    ee = cr2.r * inv(cr1.r) - eye(3);
    eta = 0.5 * (ee + ee' + ee'*ee);
    S = 1/3 * sum(eig(eta).^2);

    norms = norm(cr1.x*cr1.r - cr2.x*cr2.r,"rows");

    ## calculate dmax
    dmax = max(norms);

    ## calculate dav
    dav = mean(norms);

    ## calculate delta
    global atdb
    if (!exist("atdb","var") || isempty(atdb))
      mol_dbstart();
    endif

    deltad = (cr1.a(2)/cr1.a(1))*(cr1.a(3)/cr1.a(1))/(cr2.a(2)/cr2.a(1))/(cr2.a(3)/cr2.a(1));
    atmass = atdb.mass(cr1.ztyp);
    deltac = atmass(cr1.typ) * norms.^2 / sum(atmass(cr1.typ));
    delta = (sqrt(2) * deltac + 1) * deltad - 1;
  else
    S = dmax = dav = delta = Inf;
  endif

  ## rest
  [t1 i1 pinfo] = cr_xrd(cr1);
  i1 = i1 / trapz(t1,i1);
  [t2 i2 pinfo] = cr_xrd(cr2);
  i2 = i2 / trapz(t2,i2);
  idif = abs(cumtrapz(t1,i1-i2));
  spow = trapz(t1,idif) / (t1(end)-t1(1));

endfunction
