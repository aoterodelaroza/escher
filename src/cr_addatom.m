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

function [cr natnew typnew] = cr_addatom(cr0,x,zsy)
% function [cr natnew typnew] = cr_addatom(cr0,x,zsy)
%
% cr_addatom - add a new atom to a crystal.
%
% Input:
% cr0: the input crystal.
% x: fractional coordinates for the new atom
% zsy: either an integer (atomic number) or a string (atomic symbol/label). 
%      These are compared to cr.ztyp and cr.atyp respectively to see if 
%      the new atom falls into any of the known types
%
% Output:
% cr: the crystal containing the new atom.
% natnew: identifier for the new atom.
% typnew: type id for the new atom.
%

  cr = cr0;
  
  ## see if we have this type/create a new one
  if (ischar(zsy))
    %% do we have this atomic symbol?
    [ok idx] = ismember(zsy,cr.attyp);
    if (ok)
      nt = idx;
    else
      nt = ++cr.ntyp;
      cr.attyp(nt) = zsy;
      cr.ztyp(nt) = mol_dbatom(zsy);
    endif
  else
    %% do we have this atomic number?
    idx = find(cr.ztyp == zsy);
    if (isempty(idx))
      nt = ++cr.ntyp;
      cr.ztyp(nt) = zsy;
      cr.attyp(nt) = mol_dbsymbol(zsy);
    else
      nt = idx;
    endif
  endif
  
  ## add a new atom
  n = ++cr.nat;
  cr.typ(n) = nt;
  if (rows(3) == 3)
    x = x';
  endif
  cr.x(n,:) = x;

endfunction
