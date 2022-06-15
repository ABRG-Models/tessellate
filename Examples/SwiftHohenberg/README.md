
This directory contains code for solving the Swift-Hohenberg equations with various nonlinearities built
around the same linear operator core that ensures that the strongest growing spatially inhomogeneous
solution has a definite wavelength, in the limit a single wavelength.

It contains the following header files

region.h for creating tesselated regions
analysis.h for methods that assist with analysis of the solutions
hexGeometry.h which creates geometric objects that are compatible with hexagonal grids.


It has two solver classes
shSolver.h that solves the Swift-Hohenberg equations with real-valued fields
shPSolver.h which solves the SH equations for complex-valued fields.

these are driven by main programs
shReal.cpp for real fields
shComplex.cpp for complex fields noflux boundary conditions
shPeriodic.cpp for complex fields with periodic boundary condition

To compile mkdir build; cd build
then
     cmake ..
and
     make shReal
     make shComplex
     make shPeriodic

to make the json files

    sh shRealjson.sh
    sh shComplexjson.sh (both shComplex and shPeriodic use the same json file)

and to run

    python runshReal.py
    python runshComplex.py
    python runshPeriodic.py

NOTE: 15/06/2022 There are two different methods in shPSolver.h for identifying the contours
of zero values for a field of real values on the hexGrid, isCountourZeroPeriodic and
isContourZeroSpiral. To toggle between which one is used toggle lContourPeriodic in
shComplexjson.sh. I have included an hdf5 file first.h5 which will allow a warm start so the
solution does not need to spun up to see the effects. To use this create a directory
logsSwiftHohenberg at the level in which you run runshPeriodic.py and move first.h5 into this
directory.
