# tessellate
John Brooke's Tessellate project

This repository contains the classes for solving partial differential equations on a tesselation.

The classes are in the RegionClasses directory.

region.h creates a Voronoi tesselation based on the input of a file of seed points that is
created by a setCentres program and is called centres.inp

ksRegion.h creates tessellations based on triangles, equilateral, isosceles and scalene.

ksSolver.h is a solver class for the Keller-Segel equations.

hexGeometry.h is a class for creating geometric structures needed by the other classes.

analysis.h contains routines for the analysis of the results of solving PDEs on a tessellation.

In the Examples directory are subdirectories containing main programs for solving different
problems on tessellations.

The directories in Examples are as follows. Each directory has instructions for how to compile and
run the code for each example in a README.md file for that example. Several of the README files make
reference to the following paper, referred to in the directories as Paper 1:

@article{Brookeetal2022,
title = {Biological action at a distance: Correlated pattern formation in adjacent tessellation domains without communication},
author = {Brooke, John M. and James, Sebastian S. and Jimenez-Rodriguez, Alejandro and Wilson, Stuart P.},
journal = {PLoS computational biology},
volume = {18},
number = {3},
pages = {e1009963--e1009963}
year = {2022},
}

Examples

CookieCutter: This refers to figure 4 of Paper 1. It compares the correlations in solutions of the Keller-Segel equations on
tessellations where the boundaries do or do not have a causal influence on pattern formation.

KSInTree: This solves the Keller-Segel equations in domains which are tesselations of triangles.

MorphVis: This solves the Keller-Segel equations in domains which are Voronoi tessellations.

Schanakenberg: This solves the Schnakeberg equations in domains which are Voronoi tessllations

The remaining examples directory makes reference to the following paper


@article{Wolf2005a,
  author = {Wolf, F},
  title = {Symmetry, Multistability, and Long-Range Interactions in Brain Development},
  journal = {Physical Review Letters},
  year = {2005},
  volume = {95},
  pages = {208701-208704},
}

SwiftHohenberg: This solves the complex Swift-Hohenberg system as defined in Wolf2005a on different shaped domains.
There exists the option to set periodic or no-flux boundary conditions.

There is a also a directory Videos which contains videos of the results of the code. It currently contains
a subdirectory Subbarrels which contains videos based on the code in MorphVis and which is documented a paper currently

@article{Wolf2005a,
  author = {BrookeWilson23},
  title = {Cortical size alone cannot explain within-domain patterning},
  journal = {PLos computational biology},
  year = {2023},
  volume = {},
  pages = {},
}

this paper is based on Chapter 3 of J.M.Brooke's PhD Thesis "The effects of geometry and dynamics on biological pattern formation


being drafted
