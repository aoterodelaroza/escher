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

function [t ih pinfo] = cr_xrd(cr,root="",lambda=1.5406,fpol=0,th2limits=[5 90],npts=10001,LOG=0)
% function [t ih pinfo] = cr_xrd(cr,root="",lambda=1.5406,fpol=0,th2limits=[5 90],npts=10001,LOG=0)
%
% Powder x-ray diffractogram of a crystal (cr), for specified beam 
% wavelentgh and polarization. The output is used by xrd_gnuplot to
% do the plot.
%
% Required input variables:
% {cr}: the crystal description. The following fields are used:
%       a, b, x, typ, ztyp, qtyp.
%
% Optional input variables (all have default values):
% {lambda}: x-ray wavelength, in angstrom. The default (1.5406)
%           correspnods to the K-line of Cu.
% {fpol}: polarization of the incident line. fpol=0 (default) is 
%         unpolarized light. Use fpol=0.95 for synchrotron radiation.
% {th2limits}: 2theta-limits of the plot.
% {npoints}: number of points along the 2theta axis.
% {root}: root of the gnuplot script and data. If present a gnuplot
%         script is written. 
% {LOG}: verbose flag.
%
% Output:
% t(1:npoints): vector of 2theta values.
% ih(1:npoints): vector of intensities.
% pinfo: additional x-ray diagram information. This structure contains the 
% following fields:
% * peaks: number of detected peaks
% * th2: 2theta of the detected peaks.
% * ip: intensity of the detected peaks.
% * hvec: (hkl) label of the peaks.
% * th2limits: limits of the plot (same as in the input).
%  
  a = cr.a;
  b = cr.b;
  x = cr.x;
  typ = cr.typ;
  z = cr.ztyp;
  q = cr.qtyp;

  ieps = 1d-5;
  theps = 1d-5;
  sigma = 0.01;

  pprofile = "gaussian";
  sigma2 = sigma * sigma;
  bohr2ang = 0.52917720859;
  lambda /= bohr2ang;

  ## load atomic scattering functions
  for i = 1:length(z)
    [ascatt{i} bscatt{i} cscatt{i} c2scatt{i} Biso{i}] = dbxrd(z(i),q(i));
  endfor

  ## metric tensor et al.
  nat = size(x,1);
  cb = cos(b);
  G = [a(1)*a(1),       a(1)*a(2)*cb(3), a(1)*a(3)*cb(2);
       a(1)*a(2)*cb(3), a(2)*a(2),       a(2)*a(3)*cb(1);
       a(1)*a(3)*cb(2), a(2)*a(3)*cb(1), a(3)*a(3) ];
  omega = sqrt(det(G));
  Gr = inv(G);
  ar = sqrt(diag(Gr))';

  ## initialize intensity grid
  t = linspace(th2limits(1),th2limits(2),npts);
  ih = zeros(1,length(t));
  th2limits *= pi / 180;
  
  peaks=0;
  th2p = zeros(100,1);
  ip = zeros(100,1);
  mult = zeros(100,1);
  hvecp = zeros(100,3);
  smax = sin(th2limits(2)/2);
  hmax = ceil(2*smax/lambda/min(ar));
  for hcell = 1:hmax
    for h = -hcell:hcell
      for k = -hcell:hcell
        for l = -hcell:hcell
          if (abs(h) != hcell && abs(k) != hcell && abs(l) != hcell)
            continue
          endif
          hvec = [h k l];
          dh = sqrt(hvec * Gr * hvec');
          dh2 = dh * dh;
          dh3 = dh2 * dh;
          if (abs(lambda * dh / 2) > smax)
            continue
          endif
          th = asin(lambda * dh / 2);
          th2 = 2 * th;
          if (th2 < th2limits(1) || th > th2limits(2))
            continue
          endif
          sthlam = dh / bohr2ang / 2;
          kvec = 2 * pi * hvec;
 
          cterm = 0; sterm = 0;
          for i = 1:nat
            dh;
            if (dh < 2)
              ffac = ascatt{typ(i)}(1)*exp(-bscatt{typ(i)}(1)*dh2)+\
                     ascatt{typ(i)}(2)*exp(-bscatt{typ(i)}(2)*dh2)+\
                     ascatt{typ(i)}(3)*exp(-bscatt{typ(i)}(3)*dh2)+\
                     ascatt{typ(i)}(4)*exp(-bscatt{typ(i)}(4)*dh2)+cscatt{typ(i)};
            else
              if (z(i) == 1)
                ffac = 0;
              else
                ffac = exp(c2scatt{typ(i)}(1)+\
                           c2scatt{typ(i)}(2)*dh+\
                           c2scatt{typ(i)}(3)*dh2+\
                           c2scatt{typ(i)}(4)*dh3);
              endif
            endif
            ffac = ffac * exp(- Biso{typ(i)} * sthlam^2);
 
            cterm += ffac * cos(kvec * x(i,:)');
            sterm += ffac * sin(kvec * x(i,:)');
          endfor
          int = cterm^2 + sterm^2;
        
          ## W. Yinghua J. Appl. Cryst. 20 (1987) 258-259.
          ## profile correction
          ## lpf = (1 + cos(th2)^2) / sin(th)^2;
          ## int = int * lpf;
         
          ## gdis Lorentz-polarization correction
          ## checked in diamond
          ## same as criticized by W. Yinghua because it is the correction
          ## for the integrated intensity
          ## lpf = 0.5*(1+cos(th2)^2) / sin(th2) / sin(th);
          ## int = int * lpf;
 
          ## FoX-compatible
          ## lorentz 
          mcorr = 1 / sin(th2);
          ## slit aperture
          mcorr *= 1 / sin(th);
          ## polarization
          ## unpolarized -> fpol = 0
          ## synchrotron -> fpol = 0.95
          A = (1 - fpol) / (1 + fpol);
          mcorr *= (1 + A * (0.5 + 0.5 * cos(2*th2))) / (1+A);
          int = int * mcorr;
         
          if (int > ieps)
            if (LOG)
              printf("2theta = %.15f   [%d %d %d]  I = %.15f\n",th2*180/pi,hvec,int)
            endif
            if (strcmp(pprofile,"gaussian"))
              ih += int * exp(-(t-th2*180/pi).^2 / 2 / sigma2);
              if (all(abs(th2p-th2) > theps))
                peaks++;
                th2p(peaks) = th2;
                ip(peaks) = int;
                mult(peaks) = 1;
                hvecp(peaks,:) = hvec;
              else
                i = find(abs(th2p-th2) <= theps);
                if (length(i) > 1) 
                  [w j] = min(abs(th2p(i)-th2));
                  i = i(j);
                endif
                ## Usually the hvec with the most positive indices is the last one
                hvecp(i,:) = hvec;  
                mult(i)++;
              endif
            endif
          endif
          
        endfor
      endfor
    endfor
  endfor
  ip(1:peaks) = ip(1:peaks) .* mult(1:peaks);

  ipmax = max(ip(1:peaks));
  ihmax = max(ih);
  ih = ih / ihmax * 100;
  ip = ip / ipmax * 100;

  [th2p i j] = unique(th2p);
  ## remove the 0
  if (th2p(1) < th2limits(1))
    th2p = th2p(2:end);
    i = i(2:end);
    j -= 1;
  endif

  pinfo.peaks = length(i);
  pinfo.th2 = zeros(length(i),1);
  pinfo.ip = zeros(length(i),1);
  pinfo.hvec = zeros(length(i),3);
  for k = 1:length(i)
    pinfo.th2(k) = th2p(k) * 180 / pi;
    pinfo.ip(k) = ip(i(k));
    pinfo.hvec(k,:) = hvecp(i(k),:);
  endfor
  pinfo.th2limits = th2limits * 180 / pi;

  if (!isempty(root))
    xrd_gnuplot(root,t,ih,pinfo);
  endif

endfunction

function xrd_gnuplot(root,t,ih,pinfo)
% function xrd_gnuplot(root,t,ih,pinfo)
%
% Plot a powder x-ray diffractogram obtained using cr_xrd. This routine
% generates a data file (root.dat) and the corresponding gnuplot script
% (root.gnu).
%
% Required input variables:
% root: the root for the data file (root.dat) and gnuplot script (root.gnu).
% t: 2theta values (see cr_xrd output)
% ih: intensity values (see cr_xrd output)
% pinfo: additional peak information (see cr_xrd output).
%
% Optional input variables (all have default values):
% 
% Side effects:
% The root.dat and root.gnu files are created.
%

  # Write the data file
  [lu,msg] = fopen(strcat(root,".dat"),"w");
  if (lu == -1)
    error(sprintf("Could not open file %s.dat: %s",root,msg))
  endif
  fprintf(lu,"# 2*theta   Intensity\n")
  for i = 1:length(t)
    fprintf(lu,"%.15f %.15f\n",t(i),ih(i))
  endfor
  fclose(lu);

  # Write the gnuplot file
  [lu,msg] = fopen(strcat(root,".gnu"),"w");
  if (lu == -1)
    error(sprintf("Could not open file %s.dat: %s",root,msg))
  endif
  fprintf(lu,"set terminal postscript eps color enhanced 'Helvetica' 25\n");
  fprintf(lu,"set output '%s.eps'\n",root);
  fprintf(lu,"\n");
  for i = 1:length(pinfo.th2)
    if (pinfo.th2(i) <= pinfo.th2limits(2))
      fprintf(lu,
              "set label %d '%d%d%d' at %.6f,%.6f center rotate by 90 font 'Helvetica,12'\n",
              i, pinfo.hvec(i,:),pinfo.th2(i), min(pinfo.ip(i)+4,102));
    endif
  endfor
  fprintf(lu,"\n");
  fprintf(lu,"set xlabel '2{/Symbol Q} (degrees)'\n");
  fprintf(lu,"set ylabel 'Intensity (arb. units)'\n");
  fprintf(lu,"set xrange [%f:%f]\n",pinfo.th2limits(1),pinfo.th2limits(2));
  fprintf(lu,"set style data lines\n");
  fprintf(lu,"set grid\n");
  fprintf(lu,"unset key\n");
  fprintf(lu,"plot '%s.dat' w lines\n",root);
  fprintf(lu,"\n");
  fprintf(lu,"!epstopdf %s.eps\n",root);
  fprintf(lu,"!pdfcrop %s.pdf\n",root);
  fprintf(lu,"!mv %s-crop.pdf %s.pdf\n",root,root);
  fprintf(lu,"!rm %s.eps\n",root);
  fclose(lu);

endfunction

function [aa bb cc cc2 Biso] = dbxrd(z,q=0)
% function [aa bb cc cc2 Biso] = dbxrd(z)
%
% Coefficients for an analytical approximation to the atomic scattering factors. 
% From ITC, vol. C, 3rd ed. (2004), tables 6.1.1.4 and 6.1.1.5 (p. 578 and ff.).
%
% Required input variables:
% z: atomic number.
%
% Optional input variables:
% q: atomic charge. If a suitable q is not found in the table, defaults to zero
%    and gives a warning message.
%
% Output variables:
% {aa(4), bb(4), cc}: coefficients for the scattering factor approximation in the
% range 0 < sin(theta)/lambda < 2.0 ang^-1. The expression is:
%    f(x=sin(theta)/lambda) = sum_{i=1->4} aa(i) * exp(-bb(i) * x^2) + c
% {cc2(4)}: for the 2.0 < sin(theta)/lambda < 6.0 ang^-1, the above equation is not
% accurate, and the following is used instead:
%    ln(f(x=sin(theta)/lambda)) = sum{i=1->4} cc2(i) * x^(i-1)
% {Biso}: isotropic temperature factor.
%

  Biso = -1;
  ## H
  if (z == 1) 
    if (q == 0)
      aa = [0.489918 0.262003 0.196767 0.049879];
      bb = [20.6593 7.74039 49.5519 2.20159];
      cc = 0.001305;
      cc2 = [0. 0. 0. 0.];
      Biso = 1;
    endif
  ## C
  elseif (z == 6) 
    if (q == 0)
      aa = [2.31 1.02 1.5886 0.865];
      bb = [20.8439 10.2075 0.5687 51.6512];
      cc = 0.2156;
      cc2 = [1.70560 -1.56760 1.18930 -0.42715];
      Biso = 1;
    endif
  ## N
  elseif (z == 7) 
    if (q == 0)
      aa = [12.2126 3.13220 2.01250 1.16630];
      bb = [0.0057 9.8933 28.9975 0.5826];
      cc = -11.529;
      cc2 = [1.54940 -1.20190 0.51064 0.02472];
      Biso = 1;
    endif
  ## O
  elseif (z == 8) 
    if (q == 0)
      aa = [3.04850 2.28680 1.54630 0.867];
      bb = [13.2771 5.70110 0.3239 32.9089];
      cc = 0.2508;
      cc2 = [1.30530 -0.83742 -0.16738 0.475];
      Biso = 1;
    endif
  ## Na
  elseif (z == 11)
    if (q == 0)
      aa = [4.76260 3.17360 1.26740 1.11280];
      bb = [3.28500 8.84220 0.313600 129.424];
      cc = 0.676000;
      cc2 = [0.84558 -0.26294 -0.87884 0.76974];
      Biso = 1;
    elseif (q == 1)
      aa = [3.25650 3.93620 1.39980 1.00320];
      bb = [2.66710 6.11530 0.200100 14.0390];
      cc = 0.404000;
      cc2 = [0.84558 -0.26294 -0.87884 0.76974];
      Biso = 1;
    endif
  ## Cl
  elseif (z == 17) 
    if (q == -1)
      aa = [18.2915 7.20840 6.53370 2.33860];
      bb = [0.0066 1.1717 19.5424 60.4486];
      cc = -16.378;
      cc2 = [1.42320 -0.63936 0.84722 -0.76135];
      Biso = 1;
    elseif (q == 0)
      aa = [11.4604 7.19640 6.25560 1.64550];
      bb = [0.0104 1.16620 18.5194 47.7784];
      cc = -9.5574;
      cc2 = [1.42320 -0.63936 0.84722 -0.76135];
      Biso = 1;
    endif
  ## Co
  elseif (z == 27) 
    if (q == 0) 
      aa = [12.2841 7.34090 4.00340 2.34880];
      bb = [4.27910 0.2784 13.5359 71.1692];
      cc = 1.01180;
      cc2 = [3.959 -1.9965 3.6063 -2.3705];
      Biso = 1;
    endif
  ## Cd
  elseif (z == 48) 
    if (q == 0)
      aa = [19.2214 17.6444 4.46100 1.60290];
      bb = [0.5946 6.9089 24.7008 87.4825];
      cc = 5.0694;
      cc2 = [3.0843 -0.7145 0.84482 -0.6099];
      Biso = 1;
    endif
  else
    printf("The scattering data for atomic number %d has not been coded yet.\n")
    printf("Why not going to ITC and doing it now?\n")
    error("atom not found in xrd_data.m")
  endif  
  if (Biso == -1) 
    print_warning(q,z)
    [aa bb cc cc2 Biso] = dbxrd(z,0);
  endif
  cc2(3) /= 10;
  cc2(4) /= 100;

endfunction

function print_warning(q,z)
  printf("Warning: q = %.1f not found for Z = %.1f\n",q,z);
  printf("Either it is not coded (check ITC, vol. C!!) or not available.\n");
  printf("Using the q = 0 default\n");
endfunction

