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

function rep = rep_lightmodel(repi,model="",ifac=1)
% function rep = rep_lightmodel(repi,model)
%
% rep_lightmodel - given a representation, position the lights relative
% to the location of the camera using one of the known "models".
%
% Input variables:
% repi: input representation.
% model:
%   "tessel" - one light at the camera poistion, one in the sky
%   "tessel-noshadow" - same as tessel but lights do not cast shadows.
%   "direct" - one light at the camera, one to the side and above,
%              one to the other side and above, one behind
%   "direct-noshadow" - same as direct but no shadows.
%   "3point" - one at the side, one at the other side, one behind.
%   "3point-noshadow" - same as 3point but no shadows.
% ifac: light intensity factor.
%
% Output variables:
% rep: output representation.

  rep = repi;

  if (isnull(model))
    error("no light model given!")
  endif

  if (strcmpi(model,"tessel") || strcmpi(model,"tessel-noshadow"))
    if (strcmpi(model,"tessel"))
      ## one at the camera
      rep = rep_addlight(rep,rep.cam.location,0,ifac);
      ## one in the sky
      rep = rep_addlight(rep,rep.cam.sky*10,0,ifac);
    else
      rep = rep_addlight(rep,rep.cam.location,1,ifac);
      rep = rep_addlight(rep,rep.cam.sky*10,1,ifac);
    endif
  elseif (strcmpi(model,"direct") || strcmpi(model,"direct-noshadow"))
    ## one at the camera
    if (strcmpi(model,"direct"))
      rep = rep_addlight(rep,rep.cam.location,0,0.5*ifac);
    else
      rep = rep_addlight(rep,rep.cam.location,1,0.5*ifac);
    endif

    x0 = rep.cam.location(:);
    xc = rep.cam.lookat(:) - x0;
    dc = norm(xc); xc = xc / dc;
    xu = rep.cam.sky(:); xu = xu / norm(xu);
    xd = cross(xc,xu);
    if (norm(xd) < 1e-14)
      xd = cross(xc,rep.cam.right(:));
    endif
    xd = xd / norm(xd);

    ## fill light, 20 degrees to the side and above the camera
    thk = 20 * pi / 180;
    x = x0 + dc * tan(thk) * xd + dc * tan(thk) * xu;
    rep = rep_addlight(rep,x,1,ifac*0.3);

    ## the fill light, 20 degrees to the other side and 10 above
    ## half brightness, no shadows
    thf = 10 * pi / 180;
    x = x0 - dc * tan(thk) * xd + dc * tan(thf) * xu;
    rep = rep_addlight(rep,x,1,ifac*0.3);

    ## the rim light, opposite to the camera, no shadows
    x = x0 + 2 * dc * xc - 2 * (xc * xu') * dc * xu;
    rep = rep_addlight(rep,x,1,0.5*ifac);

  elseif (strcmpi(model,"3point") || strcmpi(model,"3point-noshadow"))
    x0 = rep.cam.location(:);
    xc = rep.cam.lookat(:) - x0;
    dc = norm(xc); xc = xc / dc;
    xu = rep.cam.sky(:); xu = xu / norm(xu);
    xd = cross(xc,xu);
    if (norm(xd) < 1e-14)
      xd = cross(xc,rep.cam.right(:));
    endif
    xd = xd / norm(xd);

    ## the key light, 20 degrees to the side and above the camera
    thk = 20 * pi / 180;
    x = x0 + dc * tan(thk) * xd + dc * tan(thk) * xu;
    if (strcmpi(model,"3point"))
      rep = rep_addlight(rep,x,0,ifac*0.5);
      rep = rep_addlight(rep,x,0,ifac*0.5);
    else
      rep = rep_addlight(rep,x,1,ifac*0.5);
      rep = rep_addlight(rep,x,1,ifac*0.5);
    endif

    ## the fill light, 20 degrees to the other side and 10 above
    ## half brightness, no shadows
    thf = 10 * pi / 180;
    x = x0 - dc * tan(thk) * xd + dc * tan(thf) * xu;
    rep = rep_addlight(rep,x,1,ifac*0.3);

    ## the rim light, opposite to the camera, no shadows
    x = x0 + 2 * dc * xc - 2 * (xc * xu') * dc * xu;
    rep = rep_addlight(rep,x,1,ifac);

  else
    error("unknown model!")
  endif

endfunction
