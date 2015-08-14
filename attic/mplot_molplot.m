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

function mplot_molplot (mol, LOG=1)
% function mplot_molplot (mol, LOG=1)
%
% mplot_molplot - add mol to the graphics database.
%
% Authors: VLC Victor Lua~na .......... <victor@carbono.quimica.uniovi.es>
%          AOR Alberto Otero-de-la-Roza <alberto@carbono.quimica.uniovi.es>
% Created: December 2011

global ge 

if (isempty(mol.conn))
   error('mol_molplot: empty connectivity matrix!');
endif

nat = length(mol.atname)
for i = 1 : nat
   ge.pt(1:3,++ge.npt) = mol.atxyz(1:3,i);
   ge.ball(++ge.nball).ct = ge.npt;
   [z,atom] = mol_dbatom(mol.atname(i),0);
   ge.ball(ge.nball).rad = atom.rad;
   ge.ball(ge.nball).rgb = atom.number;
   ge.ball(ge.nball).txt = mol.atname(i);
   point(i) = ge.npt;
   at{i} = atom;
endfor

for i = 1 : nat
   for j = i+1 : nat
      if (mol.conn(i,j))
         ge.nstk++;
         ge.stk(ge.nstk).pt1 = point(i);
         ge.stk(ge.nstk).pt2 = point(j);
         ge.stk(ge.nstk).rad = ge.defstkrad;
         if (ge.bicolor)
            ge.npt++;
            ge.pt(1:3,ge.npt) = (ge.pt(1:3,point(i))+ge.pt(1:3,point(j)))/2;
            ge.stk(ge.nstk).pt2 = ge.npt;
            ge.stk(ge.nstk).color = at{i}.number;
            ge.nstk++;
            ge.stk(ge.nstk).pt1 = ge.npt;
            ge.stk(ge.nstk).pt2 = point(j);
            ge.stk(ge.nstk).color = at{j}.number;
         else
            ge.stk(ge.nstk).color = ge.defstkcolor;
         endif
      endif
   endfor
endfor

endfunction
