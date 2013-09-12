#! /usr/bin/octave -q

addpath("../../../src");

function x = molware_sph2cart(u)
  x = [u(:,1) .* sin(u(:,2)) .* cos(u(:,3)),\
       u(:,1) .* sin(u(:,2)) .* sin(u(:,3)),\
       u(:,1) .* cos(u(:,2))];
endfunction

tic
function z = fsurf(u)
  z = 1.2*sin(pi*u(:,1)).*cos(pi*u(:,2)).*exp(-0.6*(u(:,1).*u(:,1)+u(:,2).*u(:,2)));
endfunction
function x = ftes01(u)
  x = [u(:,1),\
       u(:,2),\
       fsurf(u)];
endfunction
printf("Surface 1, tessel tes01\n");
rep = rep_surface(:,@ftes01,[-2 -2],[2 2],[101 101],:,:);
rep_write_off(rep,"surf1_none.off");
rep = rep_surface(:,@ftes01,[-2 -2],[2 2],[101 101],:,@(u)(colormode_binary(u,@fsurf)));
rep_write_coff(rep,"surf1_binary.off");
rep = rep_surface(:,@ftes01,[-2 -2],[2 2],[101 101],:,@(u)(colormode_gray(u,@fsurf)));
rep_write_coff(rep,"surf1_gray.off");
rep = rep_surface(:,@ftes01,[-2 -2],[2 2],[101 101],:,@(u)(colormode_full(u,@fsurf)));
rep_write_coff(rep,"surf1_full.off");
rep = rep_surface(:,@ftes01,[-2 -2],[2 2],[101 101],:,@(u)(colormode_hue(u,@(u)((fsurf(u)+1)*180))));
rep_write_coff(rep,"surf1_hue.off");
toc

tic
function x = ftes02(u)
  x = [ (1-0.2*cos(u(:,2))).*cos(u(:,1)),\
        (1-0.2*cos(u(:,2))).*sin(u(:,1)),\
        0.1*(sin(u(:,2))+u(:,1)/1.2-10)];
endfunction
printf("Surface 2, tessel tes02\n");
rep = rep_surface(:,@ftes02,[0, 0],[10*pi, 2*pi],[101 91],[0 1],"red");
rep_write_off(rep,"surf2.off");
toc

tic
function x = ftes03(u)
  x = [0.223015514519096 * (3*cos(u(:,1)).^2-1),\
       u(:,1),\
       u(:,2),\
       ];
  x = molware_sph2cart(x);
endfunction
printf("Surface 3, tessel tes03\n");
rep = rep_surface(:,@ftes03,[0, 0],[pi, 2*pi],[181 181],[1 1],"red");
rep_write_off(rep,"surf3.off");
toc

tic
function x = ftes04(u)
  x = [cos(u(:,2)).^3.*cos(u(:,1)).^3,\
       sin(u(:,2)).^3.*cos(u(:,1)).^3,\
       sin(u(:,1)).^3];
endfunction
printf("Surface 4, tessel tes04\n");
rep = rep_surface(:,@ftes04,[-1.4, 0],[1.4, 2*pi],[31 31],[0 0],"red");
rep_write_off(rep,"surf4.off");
toc

if (exist("legendre_sphPlm"))
  function z = rlm(u,l,m)
    z = cos(u(:,1));
    z = sqrt(2) * legendre_sphPlm(l,m,z);
    if (m > 0)
      z = z .* cos(m*u(:,2));
    elseif (m < 0)
      z = z .* sin(-m*u(:,2));
    endif 
  endfunction
  function x = frlm(u,l,m)
    x = [abs(rlm(u,l,m)),\
         u(:,1),\
         u(:,2),\
         ];
    x = molware_sph2cart(x);
  endfunction
  for l = 0:5
    for m = 0:l
      tic 
      printf("Spherical harmonic: (%d,%d)\n",l,m);
      rep = rep_surface(:,@(u)(frlm(u,l,m)),[0, 0],[pi, 2*pi],[101 101],[1 1],@(u)(colormode_binary(u,@(u)(rlm(u,l,m)))));
      rep_write_coff(rep,sprintf("ylm_%d_%d.off",l,m));
      toc
    endfor
  endfor
else
  printf("Install octave-gsl to run the spherical harmonics test.\n")
endif
