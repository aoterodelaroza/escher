
* DONE keyboard vtk window interaction (POVExport, PNG, camera position, ...)
    s key : to solid
    w key : to wireframe
* DONE NCIPlots(colourediso)
* critical points(polyhedra) from multiwfn, ORTEP(ellipsoid), popplot (read .basin)
* surface input generation for quantum espresso
* DONE (partially) mol_readxyz.m --> .py
* set camera parameters for any input data
* DONE dictionary of atom colors and covalent radii
* bond colors and double, triple bonds
* radial distribution function g(r)
* DONE use debug variable ESCHER_DEBUG='ON'
* handle units, coordinates (cartesian, cristaline, internal)
* more parsers (ORCA, ...) database of results extracted (sqlite3)
* point and space symmetry groups
* How to manage different versions of escherpy? Seamesly to C preprocessor
  conditional compilation? Different branches? How to obtain several branches
  updated? 
  git checkout cython && git rebase master |  Need to choose a master branch. How to develop in several branches?
  #git checkout master && git merge cython |
  A Cython version with a Fortran/C/C++ library to read cube files
  faster, .... Depending on available packages (Orca, Gaussian, ...) embedd or
  not package call in escher. Is pybel available?
  VTK 5 and VTK 6.
* A molden input file reader (several types exist depending on the generator package).
