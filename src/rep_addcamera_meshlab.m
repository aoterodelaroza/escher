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

function rep = rep_addcamera_meshlab(repi,t,r,persp=1,angle=45);
% function rep = rep_addcamera_meshlab(repi,t,r);
%
% rep_addcamera_meshlab - add a camera to a graphical representation
% using the trasnlation vector (1x4) and rotation vector (16x1) provided
% by meshlab. To get this numbers, open meshlab, put the camera in the
% desired orientation, then do: windows -> copy shot. Paste the contents
% of your clipboard. They will look like:
%  <!DOCTYPE ViewState>
%  <project>
%   <VCGCamera TranslationVector="-17.141 -5.59864 -3.09982 1" LensDistortion="0 0" ViewportPx="958 440" PixelSizeMm="0.0369161 0.0369161" CenterPx="479 220" FocalMm="14.0669" RotationMatrix="-0.719101 0.54336 -0.433189 0 0.285769 0.799458 0.528397 0 0.633426 0.256179 -0.730166 0 0 0 0 1 "/>
%   <ViewSettings NearPlane="1.03109" TrackScale="0.187674" FarPlane="13.0311"/>
%   <Render Lighting="1" DoubleSideLighting="0" SelectedVert="0" ColorMode="3" SelectedFace="0" BackFaceCull="0" FancyLighting="0" DrawMode="5" TextureMode="0"/>
%  </project>
%
%  The 't' and 'r' vectors are:
%  t = TranslationVector = [-17.141 -5.59864 -3.09982 1];
%  r = RotationMatrix = [-0.719101 0.54336 -0.433189 0 0.285769 0.799458 0.528397 0 0.633426 0.256179 -0.730166 0 0 0 0 1 ];
%
% angle controls the distance to the object (orthographic, persp=0) or
% the field of view (perspective, persp=1). 
%
% Required output variables:
% rep: output representation.

  ## set the camera 
  rep = repi;
  rep.cam = camera();

  pos = -t(1:3);
  rot = [r(1:3); r(5:7); r(9:11)];

  rep.cam.location = pos; # position
  rep.cam.lookat = pos + 10 * [0 0 -1] * rot; # direction
  rep.cam.persp = persp; # perspective
  rep.cam.angle = angle; # field-of-view
  rep.cam.right = [1 0 0];
  rep.cam.up = [0 0 1];
  rep.cam.sky = [0 1 0] * rot;

endfunction
