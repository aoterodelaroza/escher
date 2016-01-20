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
%       a, b, x, typ, ztyp.
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

  ieps = 1d-5;
  theps = 1d-5;
  sigma = 0.1;

  pprofile = "gaussian";
  sigma2 = sigma * sigma;
  bohr2ang = 0.52917720859;
  lambda /= bohr2ang;

  ## load atomic scattering functions
  for i = 1:length(z)
    [ascatt{i} bscatt{i} cscatt{i} c2scatt{i} Biso{i}] = dbxrd(z(i));
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
              ffac = ascatt{typ(i)}(1)*exp(-bscatt{typ(i)}(1)*dh2)+...
                     ascatt{typ(i)}(2)*exp(-bscatt{typ(i)}(2)*dh2)+...
                     ascatt{typ(i)}(3)*exp(-bscatt{typ(i)}(3)*dh2)+...
                     ascatt{typ(i)}(4)*exp(-bscatt{typ(i)}(4)*dh2)+cscatt{typ(i)};
            else
              if (z(i) == 1)
                ffac = 0;
              else
                ffac = exp(c2scatt{typ(i)}(1)+...
                           c2scatt{typ(i)}(2)*dh+...
                           c2scatt{typ(i)}(3)*dh2+...
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
            if (strcmpi(pprofile,"gaussian"))
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

function [aa bb cc cc2 Biso] = dbxrd(z)
% function [aa bb cc cc2 Biso] = dbxrd(z)
%
% Coefficients for an analytical approximation to the atomic scattering factors. 
% From ITC, vol. C, 3rd ed. (2004), tables 6.1.1.4 and 6.1.1.5 (p. 578 and ff.).
% Only neutral atoms are used.
%
% Required input variables:
% z: atomic number.
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
  sfac1 = [
           0.489918   20.6593  0.262003   7.74039  0.196767   49.5519  0.049879   2.20159  0.001305 
           0.873400   9.10370  0.630900   3.35680  0.311200   22.9276  0.178000  0.982100  0.006400 
           1.12820   3.95460  0.750800   1.05240  0.617500   85.3905  0.465300   168.261  0.037700 
           1.59190   43.6427   1.12780   1.86230  0.539100   103.483  0.702900  0.542000  0.038500 
           2.05450   23.2185   1.33260   1.02100   1.09790   60.3498  0.706800  0.140300  -0.19320 
           2.31000   20.8439   1.02000   10.2075   1.58860  0.568700  0.865000   51.6512  0.215600 
           12.2126  0.005700   3.13220   9.89330   2.01250   28.9975   1.16630  0.582600   -11.529 
           3.04850   13.2771   2.28680   5.70110   1.54630  0.323900  0.867000   32.9089  0.250800 
           3.53920   10.2825   2.64120   4.29440   1.51700  0.261500   1.02430   26.1476  0.277600 
           3.95530   8.40420   3.11250   3.42620   1.45460  0.230600   1.12510   21.7184  0.351500 
           4.76260   3.28500   3.17360   8.84220   1.26740  0.313600   1.11280   129.424  0.676000 
           5.42040   2.82750   2.17350   79.2611   1.22690  0.380800   2.30730   7.19370  0.858400 
           6.42020   3.03870   1.90020  0.742600   1.59360   31.5472   1.96460   85.0886   1.11510 
           6.29150   2.43860   3.03530   32.3337   1.98910  0.678500   1.54100   81.6937   1.14070 
           6.43450   1.90670   4.17910   27.1570   1.78000  0.526000   1.49080   68.1645   1.11490 
           6.90530   1.46790   5.20340   22.2151   1.43790  0.253600   1.58630   56.1720  0.866900 
           11.4604  0.010400   7.19640   1.16620   6.25560   18.5194   1.64550   47.7784   -9.5574 
           7.48450  0.907200   6.77230   14.8407  0.653900   43.8983   1.64420   33.3929   1.44450 
           8.21860   12.7949   7.43980  0.774800   1.05190   213.187  0.865900   41.6841   1.42280 
           8.62660   10.4421   7.38730  0.659900   1.58990   85.7484   1.02110   178.437   1.37510 
           9.18900   9.02130   7.36790  0.572900   1.64090   136.108   1.46800   51.3531   1.33290 
           9.75950   7.85080   7.35580  0.500000   1.69910   35.6338   1.90210   116.105   1.28070 
           10.2971   6.86570   7.35110  0.438500   2.07030   26.8938   2.05710   102.478   1.21990 
           10.6406   6.10380   7.35370  0.392000   3.32400   20.2626   1.49220   98.7399   1.18320 
           11.2819   5.34090   7.35730  0.343200   3.01930   17.8674   2.24410   83.7543   1.08960 
           11.7695   4.76110   7.35730  0.307200   3.52220   15.3535   2.30450   76.8805   1.03690 
           12.2841   4.27910   7.34090  0.278400   4.00340   13.5359   2.34880   71.1692   1.01180 
           12.8376   3.87850   7.29200  0.256500   4.44380   12.1763   2.38000   66.3421    1.0341 
           13.3380   3.58280   7.16760  0.247000   5.61580   11.3966   1.67350   64.8126   1.19100 
           14.0743   3.26550   7.03180  0.233300   5.16520   10.3163   2.41000   58.7097   1.30410 
           15.2354   3.06690   6.70060  0.241200   4.35910   10.7805   2.96230   61.4135   1.71890 
           16.0816   2.85090   6.37470  0.251600   3.70680   11.4468   3.68300   54.7625   2.13130 
           16.6723   2.63450   6.07010  0.264700   3.43130   12.9479   4.27790   47.7972   2.53100 
           17.0006   2.40980   5.81960  0.272600   3.97310   15.2372   4.35430   43.8163   2.84090 
           17.1789   2.17230   5.23580   16.5796   5.63770  0.260900   3.98510   41.4328   2.95570 
           17.3555   1.93840   6.72860   16.5623   5.54930  0.226100   3.53750   39.3972   2.82500 
           17.1784   1.78880   9.64350   17.3151   5.13990  0.274800   1.52920   164.934   3.48730 
           17.5663   1.55640   9.81840   14.0988   5.42200  0.166400   2.66940   132.376   2.50640 
           17.7760   1.40290   10.2946   12.8006   5.72629  0.125599   3.26588   104.354   1.91213 
           17.8765   1.27618   10.9480   11.9160   5.41732  0.117622   3.65721   87.6627   2.06929 
           17.6142   1.18865   12.0144   11.7660   4.04183  0.204785   3.53346   69.7957   3.75591 
           3.70250  0.277200   17.2356   1.09580   12.8876   11.0040   3.74290   61.6584   4.38750 
           19.1301  0.864132   11.0948   8.14487   4.64901   21.5707   2.71263   86.8472   5.40428 
           19.2674  0.808520   12.9182   8.43467   4.86337   24.7997   1.56756   94.2928   5.37874 
           19.2957  0.751536   14.3501   8.21758   4.73425   25.8749   1.28918   98.6062   5.32800 
           19.3319  0.698655   15.5017   7.98929   5.29537   25.2052  0.605844   76.8986   5.26593 
           19.2808  0.644600   16.6885   7.47260   4.80450   24.6605   1.04630   99.8156   5.17900 
           19.2214  0.594600   17.6444   6.90890   4.46100   24.7008   1.60290   87.4825   5.06940 
           19.1624  0.547600   18.5596   6.37760   4.29480   25.8499   2.03960   92.8029   4.93910 
           19.1889   5.83030   19.1005  0.503100   4.45850   26.8909   2.46630   83.9571   4.78210 
           19.6418   5.30340   19.0455  0.460700   5.03710   27.9074   2.68270   75.2825   4.59090 
           19.9644   4.81742   19.0138  0.420885   6.14487   28.5284   2.52390   70.8403   4.35200 
           20.1472   4.34700   18.9949  0.381400   7.51380   27.7660   2.27350   66.8776   4.07120 
           20.2933   3.92820   19.0298  0.344000   8.97670   26.4659   1.99000   64.2658   3.71180 
           20.3892   3.56900   19.1062  0.310700   10.6620   24.3879   1.49530   213.904   3.33520 
           20.3361   3.21600   19.2970  0.275600   10.8880   20.2073   2.69590   167.202   2.77310 
           20.5780   2.94817   19.5990  0.244475   11.3727   18.7726   3.28719   133.124   2.14678 
           21.1671   2.81219   19.7695  0.226836   11.8513   17.6083   3.33049   127.113   1.86264 
           22.0440   2.77393   19.6697  0.222087   12.3856   16.7669   2.82428   143.644   2.05830 
           22.6845   2.66248   19.6847  0.210628   12.7740   15.8850   2.85137   137.903   1.98486 
           23.3405   2.56270   19.6095  0.202088   13.1235   15.1009   2.87516   132.721   2.02876 
           24.0042   2.47274   19.4258  0.196451   13.4396   14.3996   2.89604   128.007   2.20963 
           24.6274   2.38790   19.0886  0.194200   13.7603   13.7546   2.92270   123.174   2.57450 
           25.0709   2.25341   19.0798  0.181951   13.8518   12.9331   3.54545   101.398   2.41960 
           25.8976   2.24256   18.2185  0.196143   14.3167   12.6648   2.95354   115.362   3.58324 
           26.5070   2.18020   17.6383  0.202172   14.5596   12.1899   2.96577   111.874   4.29728 
           26.9049   2.07051   17.2940  0.197940   14.5583   11.4407   3.63837   92.6566   4.56796 
           27.6563   2.07356   16.4285  0.223545   14.9779   11.3604   2.98233   105.703   5.92046 
           28.1819   2.02859   15.8851  0.238849   15.1542   10.9975   2.98706   102.961   6.75621 
           28.6641   1.98890   15.4345  0.257119   15.3087   10.6647   2.98963   100.417   7.56672 
           28.9476   1.90182   15.2208   9.98519   15.1000  0.261033   3.71601   84.3298   7.97628 
           29.1440   1.83262   15.1726   9.59990   14.7586  0.275116   4.30013   72.0290   8.58154 
           29.2024   1.77333   15.2293   9.37046   14.5135  0.295977   4.76492   63.3644   9.24354 
           29.0818   1.72029   15.4300   9.22590   14.4327  0.321703   5.11982   57.0560   9.88750 
           28.7621   1.67191   15.7189   9.09227   14.5564  0.350500   5.44174   52.0861   10.4720 
           28.1894   1.62903   16.1550   8.97948   14.9305  0.382661   5.67589   48.1647   11.0005 
           27.3049   1.59279   16.7296   8.86553   15.6115  0.417916   5.83377   45.0011   11.4722 
           27.0059   1.51293   17.7639   8.81174   15.7131  0.424593   5.78370   38.6103   11.6883 
           16.8819  0.461100   18.5913   8.62160   25.5582   1.48260   5.86000   36.3956   12.0658 
           20.6809  0.545000   19.0417   8.44840   21.6575   1.57290   5.96760   38.3246   12.6089 
           27.5446  0.655150   19.1584   8.70751   15.5380   1.96347   5.52593   45.8149   13.1746 
           31.0617  0.690200   13.0637   2.35760   18.4420   8.61800   5.96960   47.2579   13.4118 
           33.3689  0.704000   12.9510   2.92380   16.5877   8.79370   6.46920   48.0093   13.5782 
           34.6726  0.700999   15.4733   3.55078   13.1138   9.55642   7.02588   47.0045   13.6770 
           35.3163  0.685870   19.0211   3.97458   9.49887   11.3824   7.42518   45.4715   13.7108 
           35.5631  0.663100   21.2816   4.06910   8.00370   14.0422   7.44330   44.2473   13.6905 
           35.9299  0.646453   23.0547   4.17619   12.1439   23.1052   2.11253   150.645   13.7247 
           35.7630  0.616341   22.9064   3.87135   12.4739   19.9887   3.21097   142.325   13.6211 
           35.6597  0.589092   23.1032   3.65155   12.5977   18.5990   4.08655   117.020   13.5266 
           35.5645  0.563359   23.4219   3.46204   12.7473   17.8309   4.80703   99.1722   13.4314 
           35.8847  0.547751   23.2948   3.41519   14.1891   16.9235   4.17287   105.251   13.4287 
           36.0228  0.529300   23.4128   3.32530   14.9491   16.0927   4.18800   100.613   13.3966 
           36.1874  0.511929   23.5964   3.25396   15.6402   15.3622   4.18550   97.4908   13.3573 
           36.5254  0.499384   23.8083   3.26371   16.7707   14.9455   3.47947   105.980   13.3812 
           36.6706  0.483629   24.0992   3.20647   17.3415   14.3136   3.49331   102.273   13.3592 
           36.6488  0.465154   24.4096   3.08997   17.3990   13.4346   4.21665   88.4834   13.2887 
           36.7881  0.451018   24.7736   3.04619   17.8919   12.8946   4.23284   86.0030   13.2754 
           36.9185  0.437533   25.1995   3.00775   18.3317   12.4044   4.24391   83.7881   13.2674 
           ];
  sfac2 = [
           0.0  0.0            0.0  0.0       1.0000 
           0.52543  -3.43300   4.80070  -2.54760  1.0000 
           0.89463  -2.43660   2.32500  -0.71949  1.0000 
           1.25840  -1.94590   1.30460  -0.04297  1.0000 
           1.66720  -1.85560   1.60440  -0.65981  1.0000 
           1.70560  -1.56760   1.18930  -0.42715  1.0000 
           1.54940  -1.20190   0.51064  0.02472   1.0000 
           1.30530  -0.83742  -0.16738  0.47500   1.0000 
           1.16710  -0.63203  -0.40207  0.54352   1.0000 
           1.09310  -0.50221  -0.53648  0.60957   0.9995 
           0.84558  -0.26294  -0.87884  0.76974   1.0000 
           0.71877  -0.13144  -1.20900  0.82738   1.0000 
           0.67975  -0.08756  -0.95431  0.72294   1.0000 
           0.70683  -0.09888  -0.98356  0.55631   1.0000 
           0.85532  -0.21262  -0.37390  0.20731   1.0000 
           1.10400  -0.40325   0.20094  -0.26058  1.0000 
           1.42320  -0.63936   0.84722  -0.76135  0.9995 
           1.82020  -0.92776   1.59220  -1.32510  0.9995 
           2.26550  -1.24530   2.38330  -1.91290  0.9990 
           2.71740  -1.55670   3.13170  -2.45670  0.9990 
           3.11730  -1.81380   3.71390  -2.85330  0.9990 
           3.45360  -2.01150   4.13170  -3.11710  0.9995 
           3.71270  -2.13920   4.34610  -3.22040  0.9995 
           3.87870  -2.19000   4.38670  -3.17520  1.0000 
           3.98550  -2.18850   4.27960  -3.02150  1.0000 
           3.99790  -2.11080   3.98170  -2.71990  1.0000 
           3.95900  -1.99650   3.60630  -2.37050  1.0000 
           3.86070  -1.88690   3.12390  -1.94290  1.0000 
           3.72510  -1.65500   2.60290  -1.49760  0.9995 
           3.55950  -1.45100   2.03390  -1.02160  0.9995 
           3.37560  -1.23910   1.46160  -0.55471  0.9995 
           3.17800  -1.02230   0.89119  -0.09884  0.9995 
           2.97740  -0.81038   0.34861  0.32231   0.9995 
           2.78340  -0.61110  -0.14731  0.69837   0.9995 
           2.60610  -0.43308  -0.57381  1.00950   0.9995 
           2.44280  -0.27244  -0.95570  1.27070   0.9995 
           2.30990  -0.14328  -1.22600  1.45320   1.0000 
           2.21070  -0.04770  -1.41100  1.55410   1.0000 
           2.14220  0.01935   -1.52240  1.59630   1.0000 
           2.12690  0.08618   -1.49190  1.51820   1.0000 
           2.12120  0.05381   -1.50070  1.50150   1.0000 
           2.18870  -0.00655  -1.25340  1.24010   1.0000 
           2.25730  -0.05737  -1.07450  1.06630   1.0000 
           2.37300  -0.15040  -0.77694  0.79060   0.9995 
           2.50990  -0.25906  -0.44719  0.49443   0.9995 
           2.67520  -0.39137  -0.05894  0.15404   0.9995 
           2.88690  -0.56119   0.42189  -0.25659  0.9990 
           3.08430  -0.71450   0.84482  -0.60990  0.9990 
           3.31400  -0.89697   1.35030  -1.03910  0.9990 
           3.49840  -1.02990   1.68990  -1.29860  0.9990 
           3.70410  -1.18270   2.08920  -1.61640  0.9990 
           3.88240  -1.30980   2.41170  -1.86420  0.9990 
           4.08010  -1.45080   2.76730  -2.13920  0.9990 
           4.24610  -1.56330   3.04200  -2.34290  0.9990 
           4.38910  -1.65420   3.25450  -2.49220  0.9995 
           4.51070  -1.72570   3.41320  -2.59590  0.9995 
           4.60250  -1.77070   3.49970  -2.64050  0.9995 
           4.69060  -1.81790   3.60280  -2.70670  0.9995 
           4.72150  -1.81390   3.56480  -2.65180  0.9995 
           4.75090  -1.80800   3.51970  -2.59010  1.0000 
           4.74070  -1.76600   3.37430  -2.44210  1.0000 
           4.71700  -1.71410   3.20800  -2.28170  1.0000 
           4.66940  -1.64140   2.98580  -2.07460  1.0000 
           4.61010  -1.55750   2.73190  -1.84040  0.9995 
           4.52550  -1.45520   2.43770  -1.57950  0.9995 
           4.45230  -1.36440   2.17540  -1.34550  0.9990 
           4.37660  -1.27460   1.92540  -1.13090  0.9990 
           4.29460  -1.18170   1.67060  -0.91467  0.9990 
           4.21330  -1.09060   1.42390  -0.70804  0.9990 
           4.13430  -1.00310   1.18810  -0.51120  0.9990 
           4.04230  -0.90518   0.92889  -0.29820  0.9990 
           3.95160  -0.80978   0.68951  -0.09620  0.9990 
           3.85000  -0.70599   0.41103  0.11842   0.9990 
           3.76510  -0.61807   0.18568  0.29787   0.9990 
           3.67600  -0.52688  -0.04706  0.48180   0.9995 
           3.60530  -0.45420  -0.22529  0.61700   0.9995 
           3.53130  -0.37856  -0.41174  0.75967   0.9995 
           3.47070  -0.31534  -0.56487  0.87492   0.9995 
           3.41630  -0.25987  -0.69030  0.96224   0.9995 
           3.37350  -0.21428  -0.79013  1.02850   1.0000 
           3.34590  -0.18322  -0.84911  1.05970   1.0000 
           3.32330  -0.15596  -0.89878  1.08380   1.0000 
           3.31880  -0.14554  -0.90198  1.06850   1.0000 
           3.32030  -0.13999  -0.89333  1.04380   1.0000 
           3.34250  -0.15317  -0.83350  0.97641   1.0000 
           3.37780  -0.17800  -0.74320  0.88510   1.0000 
           3.41990  -0.20823  -0.64000  0.78354   0.9995 
           3.47530  -0.25005  -0.50660  0.65836   0.9995 
           3.49020  -0.25109  -0.49651  0.64340   0.9995 
           3.61060  -0.35409  -0.18926  0.36849   0.9995 
           3.68630  -0.41329  -0.01192  0.20878   0.9995 
           3.76650  -0.47542   0.16850  0.05060   0.9990 
           3.82870  -0.51955   0.29804  -0.06566  0.9990 
           3.88970  -0.56296   0.42597  -0.18080  0.9990 
           3.95060  -0.60554   0.54967  -0.29112  0.9985 
           4.01470  -0.65062   0.67922  -0.40588  0.9985 
           4.07780  -0.69476   0.80547  -0.51729  0.9985 
           4.14210  -0.73977   0.93342  -0.62981  0.9980 
           ];
  if (z < 1 || z > 98) 
    error("atomic scattering factor unknown")
  endif

  aa = [sfac1(z,1) sfac1(z,3) sfac1(z,5) sfac1(z,7)];
  bb = [sfac1(z,2) sfac1(z,4) sfac1(z,6) sfac1(z,8)];
  cc = sfac1(z,9);
  cc2 = [sfac2(z,1), sfac2(z,2), sfac2(z,3)/10, sfac2(z,4)/100];
  Biso = sfac2(z,5);

endfunction

