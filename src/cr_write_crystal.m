% Copyright (C) 2012 Victor Lua~na and Alberto Otero-de-la-Roza
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

function cr_write_crystal(cr,file="",basisfile="",acpfile="")
% function cr_write_crystal(cr,file="",basisfile="",acpfile="")
%
% cr_write_crysatl -- write a basic crystal input file.
%
% Input:
% cr: crystal structure.
% file: name of the output file. If no file is given (default), then
%       the file is written to the stdout.
%

  bohrtoans = 0.52917720859;
  if (!isstruct(cr) || isempty(cr))
    error("Invalid cr input")
  endif

  cr_write_critic2(cr,file);
  [s out] = system(sprintf("critic2 %s",file));
  if (s != 0)
    return
  endif
  alines = strsplit(out,"\n");

  spg = -1;
  nat = 0;
  atpos = [];
  atz = [];
  holo = "";
  for i = 1:length(alines)-1
    alines{i};
    if (spg < 0 && regexp(alines{i},"Space group.*Hermann-Mauguin.*"))
      aa = strsplit(alines{i});
      spg = str2num(aa{end}(1:end-1));
    elseif (length(holo) == 0 && regexp(alines{i},"Holohedry"))
      holo = strsplit(alines{i});
      holo = holo{3};
    elseif (nat == 0 && regexp(alines{i},"List of non-equivalent atoms in the unit cell"))
      j = i+1;
      while (1)
        j++;
        if (regexp(alines{j},"List"))
          break
        endif
        nat++;
        [id x y z at mult iz] = sscanf(alines{j},"%d %f %f %f %s %d %d","C");
        atpos(nat,1) = x;
        atpos(nat,2) = y;
        atpos(nat,3) = z;
        atz(nat) = iz;
      endwhile
    endif
  endfor
  if (spg < 0)
    error("Space group not found.");
  endif
  if (nat == 0)
    error("No atoms found.");
  endif
  if (length(holo) == 0 || strcmpi(holo,"unknown"))
    error("No holohedry found.");
  endif

  if (!isempty(file))
    lu = fopen(file,"w");
  else
    lu = stdout();
  endif

  aa = cr.a * bohrtoans;
  bb = cr.b * 180 / pi;
  if (strcmpi(holo,"trigonal") && (any(abs(bb(1:2)-90) > 1e-8) || abs(bb(3)-120) > 1e-8))
    rhomb = 1;
  else
    rhomb = 0;
  endif

  fprintf(lu,"Title\n");
  fprintf(lu,"CRYSTAL\n");
  fprintf(lu,"0 %d 0\n",rhomb);
  fprintf(lu,"%d\n",spg);

  if (strcmpi(holo,"triclinic"))
    fprintf(lu,"%.10f %.10f %.10f %.4f %.4f %.4f\n",aa,bb);
  elseif (strcmpi(holo,"monoclinic"))
    fprintf(lu,"%.10f %.10f %.10f %.4f\n",aa,bb(find(abs(bb - 90) > 1e-8)));
  elseif (strcmpi(holo,"orthorhombic"))
    fprintf(lu,"%.10f %.10f %.10f\n",aa);
  elseif (strcmpi(holo,"tetragonal"))
    fprintf(lu,"%.10f %.10f\n",aa(1),aa(3));
  elseif (strcmpi(holo,"trigonal"))
    if (rhomb)
      fprintf(lu,"%.10f %.4f\n",aa(1),bb(1));
    else
      fprintf(lu,"%.10f %.10f\n",aa(1),aa(3));
    endif
  elseif (strcmpi(holo,"hexagonal"))
    fprintf(lu,"%.10f %.10f\n",aa(1),aa(3));
  elseif (strcmpi(holo,"cubic"))
    fprintf(lu,"%.10f\n",aa(1));
  else
    error("Unknown holohedry.");
  endif

  if (length(acpfile) > 0 && exist(acpfile,"file"))
    izadd = 200;
    acp = parseacp(acpfile);
  else
    izadd = 0;
  endif
  fprintf(lu,"%d\n",nat);
  for i = 1:nat
    fprintf(lu,"%d %.10f %.10f %.10f\n",izadd+atz(i),atpos(i,:));
  endfor

  fprintf(lu,"SETPRINT\n");
  fprintf(lu,"1\n");
  fprintf(lu,"72 1\n");
  fprintf(lu,"OPTGEOM\n");
  fprintf(lu,"FULLOPTG\n");
  fprintf(lu,"TOLDEE\n");
  fprintf(lu,"7\n");
  fprintf(lu,"TOLDEG\n");
  fprintf(lu,"0.0003\n");
  fprintf(lu,"TOLDEX\n");
  fprintf(lu,"0.0012\n");
  fprintf(lu,"ENDOPT\n");
  fprintf(lu,"END\n");

  if (!isempty(basisfile))
    basis = parsebasis(basisfile);
  else
    error("Basis file not provided.")
  endif
  izalist = unique(sort(atz));
  for i = 1:length(izalist)
    iz = izalist(i);
    sym = mol_dbsymbol(iz);
    
    found = 0;
    for j = 1:length(basis)
      if (strcmpi(basis{j}.atom,sym))
        fprintf(lu,"%d %d\n",iz+izadd,basis{j}.nblock);
        if (length(acpfile) > 0 && exist(acpfile,"file"))
          found2 = 0;
          for k = 1:length(acp)
            if (strcmpi(acp{k}.atom,sym))
              this = acp{k};
              found2 = 1;
            endif
          endfor
          if (!found2)
            error(sprintf("Atom %s not found in acp file",sym))
          endif
          fprintf(lu,"INPUT\n");

          alist = {"l","s","p","d"};
          nterm = [0 0 0 0];
          nout = 0;
          aout = [];
          lastk = 1;
          lastnout = 0;
          for k = 1:length(alist)
            found2 = 0;
            for l = 1:this.nblock
              if (strcmpi(this.block{l}.name,alist{k}))
                found2 = 1;
                nterm(k) += this.block{l}.nterm;
                for m = 1:this.block{l}.nterm
                  nout++;
                  aout(nout,1) = this.block{l}.exp(m);
                  aout(nout,2) = this.block{l}.coef(m);
                endfor
                lastk = k;
                lastnout = nout;
              endif
            endfor
            if (!found2)
              nterm(k) = 1;
              nout++;
              aout(nout,1) = 20.;
              aout(nout,2) = 0.;
            endif
          endfor

          nterm(lastk+1:4) = 0;
          fprintf(lu,"%d. %d %d %d %d 0 0\n",iz,nterm);
          for k = 1:lastnout
            fprintf(lu,"%.10f %.14f 0\n",aout(k,1),aout(k,2));
          endfor
        endif
        iztot = 0;
        for k = 1:basis{j}.nblock
          if (strcmpi(basis{j}.block{k}.type,"s"))
            ild = 0;
            ival = 2;
          elseif (strcmpi(basis{j}.block{k}.type,"p"))
            ild = 2;
            ival = 6;
          elseif (strcmpi(basis{j}.block{k}.type,"d"))
            ild = 3;
            ival = 10;
          elseif (strcmpi(basis{j}.block{k}.type,"f"))
            ild = 4;
            ival = 14;
          endif
          if (k != basis{j}.nblock)
            iztot += ival;
          else
            ival = iz - iztot;
          endif
          fprintf(lu,"0 %d %d %d. 0.\n",ild,basis{j}.block{k}.nprim,ival);
          for l = 1:basis{j}.block{k}.nprim
            fprintf(lu,"%.10f %.10f\n",basis{j}.block{k}.exp(l),...
                    basis{j}.block{k}.coef1(l));
          endfor
        endfor
        found = 1;
        break
      endif
    endfor
    if (!found)
      error(sprintf("Basis for atom %s not found",sym))
    endif
  endfor
  fprintf(lu,"99 0\n");
  fprintf(lu,"END\n");
  fprintf(lu,"GRIMME\n");
  fprintf(lu,"1.00 20. 25.\n");
  fprintf(lu,"4\n");
  fprintf(lu,"%d 0.14 1.001\n",1+izadd);
  fprintf(lu,"%d 1.75 1.452\n",6+izadd);
  fprintf(lu,"%d 1.23 1.397\n",7+izadd);
  fprintf(lu,"%d 0.70 1.342\n",8+izadd);
  fprintf(lu,"SHRINK\n");
  fprintf(lu,"4 4\n");
  fprintf(lu,"TOLDEE\n");
  fprintf(lu,"7\n");
  fprintf(lu,"END\n");

  if (!isempty(file))
    fclose(lu);
  endif

endfunction

function basis = parsebasis(alist)
  global ferr

  if (ferr > 0) 
    fprintf(ferr,"# Start parsebasis - %s\n",strtrim(ctime(time())));
    fflush(ferr);
  endif

  ## Accept a string instead of a cell array
  if (ischar(alist))
    if (!exist(alist,"file"))
      basis = alist;
      return
    else
      alist = {alist};
    endif
  endif
  
  ## Run over all files and read the basis files
  nbas = 0;
  basis = cell();
  for i = 1:length(alist)
    ## Debug
    if (ferr > 0)
      fprintf(ferr,"# Reading file %d (%s) at %s\n",i,alist{i},strtrim(ctime(time())));
      fflush(ferr);
    endif

    ## Check that the file exists
    if (!exist(alist{i},"file"))
      error(sprintf("Initial basis file not found: %s\n",alist{i}));
    endif
    ## Open and read the file
    fid = fopen(alist{i},"r");
    if (fid < 0)
      error(sprintf("Error opening file: %s",files{i}));
    endif
    while (!feof(fid))
      ## Skip blank lines and four-star lines
      line = strtrim(fgetl(fid));
      if (!ischar(line) || length(line) == 0)
        continue
      end
      if (strcmp(line,"****"))
        continue
      endif
      ## This is the 'Atom 0' line at the beginning
      curatoms = {};
      anew = strsplit(line);
      curatoms = anew(1:end-1);
      for j = 1:length(curatoms)
        curatoms{j} = strrep(curatoms{j},"-","");
      endfor
      ## Add the new basis set for this atom
      nbas++;
      basis{nbas}.atom = curatoms{1};
      basis{nbas}.name = curatoms{1};
      basis{nbas}.nblock = 0;
      basis{nbas}.block = cell();
      line = strtrim(fgetl(fid));
      if (!ischar(line) || length(line) == 0)
        error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
      end
      while (!strcmp(line,"****"))
        basis{nbas}.nblock++;
        [ityp inum rdum] = sscanf(line,"%s %d %f","C");
        if (isempty(inum))
          basis{nbas}.bstring = ityp;
          basis{nbas}.nblock = 0;
          break
        else
          basis{nbas}.bstring = [];
          basis{nbas}.block{basis{nbas}.nblock}.type = ityp;
          basis{nbas}.block{basis{nbas}.nblock}.bit1 = rdum;
          basis{nbas}.block{basis{nbas}.nblock}.nprim = inum;
          basis{nbas}.block{basis{nbas}.nblock}.exp = zeros(1,inum);
          basis{nbas}.block{basis{nbas}.nblock}.coef1 = zeros(1,inum);
          basis{nbas}.block{basis{nbas}.nblock}.coef2 = zeros(1,inum);
        endif
        for k = 1:inum
          line = strtrim(fgetl(fid));
          if (!ischar(line) || length(line) == 0)
            error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
          end
          [rexp rcoef rcoef2] = sscanf(line,"%f %f %f","C");
          basis{nbas}.block{basis{nbas}.nblock}.exp(k) = rexp;
          basis{nbas}.block{basis{nbas}.nblock}.coef1(k) = rcoef;
          if (isempty(rcoef2))
            basis{nbas}.block{basis{nbas}.nblock}.coef2(k) = NaN;
          else
            basis{nbas}.block{basis{nbas}.nblock}.coef2(k) = rcoef2;
          endif
        endfor
        line = strtrim(fgetl(fid));
        if (!ischar(line) || length(line) == 0)
          error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
        end
      endwhile

      nbas0 = nbas;
      for k = 2:length(curatoms)
        nbas++;
        basis{nbas} = basis{nbas0};
        basis{nbas}.atom = curatoms{k};
        basis{nbas}.name = curatoms{k};
      endfor
    endwhile
    ## Debug
    if (ferr > 0)
      fprintf(ferr,"# Read basis set for atoms: ");
      for i = 1:nbas
        fprintf(ferr,"%s ",basis{i}.atom);
      endfor
      fprintf(ferr,"\n");
      fflush(ferr);
    endif
    fclose(fid);
  endfor

  if (ferr > 0) 
    fprintf(ferr,"# End parsebasis - %s\n",strtrim(ctime(time())));
    fflush(ferr);
  endif

end

function acp = parseacp(alist)
  %% function acp = parseacp(alist)
  %%
  %% Read the alist file containing a ACP specification in Gaussian
  %% format for one or more atom.  alist can be either a string or a
  %% cell structure of strings. All strings should correspond to a
  %% file. The parseacp routine returns the ACPs in a cell array of
  %% structures.

  global ferr

  ## Debug
  if (ferr > 0) 
    fprintf(ferr,"# Start parseacp - %s\n",strtrim(ctime(time())));
    fflush(ferr);
  endif

  ## Accept a string instead of a cell array
  if (ischar(alist))
    alist = {alist};
  endif
  
  ## Run over all files and read the ACPs
  nacp = 0;
  acp = cell();
  for i = 1:length(alist)
    ## Check that the file exists
    if (!exist(alist{i},"file"))
      error(sprintf("Initial ACP file not found: %s\n",alist{i}));
    endif
    ## Open and read the file
    fid = fopen(alist{i},"r");
    if (fid < 0)
      error(sprintf("Error opening file: %s",files{i}));
    endif
    while (!feof(fid))
      line = fgetl(fid);
      if (feof(fid))
        break
      endif
      ## Skip blank lines and four-star lines
      line = strtrim(line);
      if (!ischar(line) || length(line) == 0)
        continue
      end
      if (strcmp(line,"****"))
        continue
      endif
      ## This is the 'Atom 0' line at the beginning. Add the new ACP
      anew = strsplit(line);
      nacp++;
      acp{nacp}.atom = strrep(anew{1},"-","");
      ## Read the second line
      line = strtrim(fgetl(fid));
      if (!ischar(line) || length(line) == 0)
        error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
      end
      anew = strsplit(line);
      if (length(anew) != 3)
        error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
      endif
      acp{nacp}.name = anew{1};
      acp{nacp}.nblock = str2num(anew{2})+1;
      acp{nacp}.block = cell(1,acp{nacp}.nblock);
      acp{nacp}.nelec = str2num(anew{3});
      acp{nacp}.nparam = 0;
      ## Read all the blocks
      for i = 1:acp{nacp}.nblock
        ## Title of the block
        line = strtrim(fgetl(fid));
        if (!ischar(line) || length(line) == 0)
          error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
        end
        acp{nacp}.block{i}.name = line;

        ## Number of terms in the ACP
        line = strtrim(fgetl(fid));
        if (!ischar(line) || length(line) == 0)
          error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
        end
        acp{nacp}.block{i}.nterm = sscanf(line,"%d","C");
        acp{nacp}.block{i}.n = zeros(acp{nacp}.block{i}.nterm,1);
        acp{nacp}.block{i}.exp = zeros(acp{nacp}.block{i}.nterm,1);
        acp{nacp}.block{i}.coef = zeros(acp{nacp}.block{i}.nterm,1);
        for j = 1:acp{nacp}.block{i}.nterm
          line = strtrim(fgetl(fid));
          if (length(line) == 0)
            error(sprintf("Incorrect Gaussian ECP format in file %s, line '%s'",alist{i},line));
          end
          [a1 a2 a3] = sscanf(line,"%f %f %f","C");
          acp{nacp}.block{i}.n(j) = a1;
          acp{nacp}.block{i}.exp(j) = a2;
          acp{nacp}.block{i}.coef(j) = a3;
          acp{nacp}.nparam += 2;
        endfor
      endfor
    endwhile
    fclose(fid);
  endfor

  ## Debug
  if (ferr > 0) 
    fprintf(ferr,"# End parseacp - %s\n",strtrim(ctime(time())));
    fflush(ferr);
  endif
end
