
This directory contains code for solving the Swift-Hohenberg equations with various nonlinearities built
around the same linear operator core that ensures that the strongest growing spatially inhomogeneous
solution has a definite wavelength, in the limit a single wavelength.

It contains the following header files

region.h for creating tesselated regions
analysis.h for methods that assist with analysis of the solutions
hexGeometry.h which creates geometric objects that are compatible with hexagonal grids.


It has two solver classes
shSolver.h that solves the Swift-Hohenberg equations with real-valued fields
shCSolver.h which solves the SH equations for complex-valued fields.

each of these has a main program
sHreal.cpp for real fields
sHcomplex.cpp for complex fields

To compile mkdir build; cd build
then
     cmake ..
and
     make sHreal
     make sHcomplex

to make the json files

    sh sHrealjson.sh
    sh sHcomplexjson.sh

and to run

    python runsHreal.py
    python runsHcomplex.py
