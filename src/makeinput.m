#! /usr/bin/octave-cli -q

## default values for output format
format long;
usestdout = 0;
fun="blyp";

## Help function and documentation for this routine
function helpme()
  helpstr = "\n\
Utility to transform input and output formats. ##\n\
\n\
syntax:\n\
  makeinput.m input.xxx output.yyy option1 option2 option3 ...\n\
\n\
where the format is detected by the extension. \n\
  input: xyz.\n\
  output: \n\
   gjf - Gaussian09\n\
   turbo - turbomole\n\
\n\
options:\n\
   - or stdout : write to the standard output.\n\
   bsic : write a turbomole define template for use with BSIC\n\
";
  disp(helpstr);
endfunction

## find the name and extension of a filename
function [name ext] = name_and_extension(str)
  idx = index(str,".","last");
  name = str(1:idx-1);
  ext = str(idx+1:end);
endfunction

## Get the list of parameters passed to makeinput
arg_list = argv();
if (nargin < 2) 
  helpme()
  return
endif

## process the input file
file = arg_list{1};
if (!exist(file,"file"))
  error(sprintf("Could not find input file: %s",file))
endif
[name ext] = name_and_extension(file);
if (strcmpi(ext,"xyz")) 
  mol = mol_readxyz(file);
else
  error(sprintf("Unknown filename extension: %s",file))
endif

## check the sanity of the molecule
if (mol.nat <= 0)
  error(sprintf("Error reading molecule in file: %s",file))
endif

## process the output file
file = arg_list{2};
[name ext] = name_and_extension(file);
if (strcmpi(ext,"gjf") || strcmpi(ext,"turbo")) 
  oformat = ext;
else
  error(sprintf("Unknown filename extension: %s",file))
endif

## process the remaining options
for i = 3:nargin
  opt = arg_list{i};
  if (strcmpi(opt,"-") || strcmpi(opt,"stdout"))
    usestdout = 1;
  elseif (strcmpi(opt,"bsic"))
    turbobsic = 1;
  elseif (length(strtrim(opt)) > 0)
    error(sprintf("Unknown option: %s",opt))
  else
    break
  endif
endfor

## open the output file
if (usestdout) 
  fid = stdout();
  file = "none";
else
  fid = fopen(file,"w");
endif

## write the file using escher's routines
err = 0;
if (strcmp(oformat,"gjf"))
  err = mol_writegjf(mol,file);
elseif (strcmp(oformat,"turbo"))
  if (!exist("turbobsic","var"))
    turbobsic = 0;
  endif
  mol_writeturbo(mol,file,turbobsic);
else
  error(sprintf("Unknown format: %s",oformat))
endif
if (err) 
  error(sprintf("Problem writing file: %s",file))
endif

## close the output file
if (!usestdout)
  fclose(fid);
endif

