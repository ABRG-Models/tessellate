
#!/bin/zsh

## Set the parameters for displayCircles, to apply the cookie cutter method
## and compare it with the solutions that evolve in polygons
## logpath : directory for output
## dt: timestep for integration
## Dn: diffusion constant for n equation
## Dchi: chemotaxis constant
## Dc: diffusion constant for c equation
## scale: controls the hex-to-hex resolution of the grid
## xspan: controls the width of the grid
## numsteps: maximum number of integration steps
## numAdjust: interval to reset boundary conditions if > numsteps is ignored
## numprint: not used in this code
## Lcontinue: if false start from random i.c., else start from previous .h5 files
## LfixedSeed: if false take random seed from clock, if true use a value fixed in the code
## Lgraphics: if true then output graphics, false is used when running remotely.
## LDn: if true then Dn is adjusted to the size of each region, if false it is uniform
## numSectors: number of sectors used in evaluating the angular and radial degrees
## aNoiseGain: controls the strength of noise
## boundaryFalloffDist: defines boundary layer when creating random i.c.
## nnInitialOffset: value of n for spatially uniform solution
## ccInitialOffset: value of c for spatially uniform solution
## overwrite_logs: not used in this code
## skipMorph: if true return after the morph 0 run, i.e. just calculate for polygonal Voronoi tessellation
## lPerturb: not used in this code
## iter: not used in this code
## lminRadius: not used in this code
## off: not used in this code
## plotevery: controls frequency of plotting
## checkevery: controls frequency of convergence checks
## saveplots: if true save plots as .png files
## vidframes: if true set plotnames for video creation
## win_width: controls size of plotting window
## fov how much we shrink or magnify the view
## NUMPOINTS: number of seed points for Voronoi tessellation
## diffTol tolerance for convergence
## lengthScale the lengthscale for the Laplacian operators


            cat > ./pFieldVis.json <<EOF
{
    "logpath" : "./logsMorph",
    "dt": 0.00005,
    "Dn": $1,
    "Dchi" : $2,
    "Dc" : $3,
    "scale": 8,
    "xspan": 5.0,
    "numsteps": 300000,
    "numprint": 999,
    "numAdjust" : 10000,
    "Lcontinue" : false,
    "LfixedSeed" : false,
    "Lgraphics" : true,
    "LDn" : false,
    "lBoundZero" : true,
    "numSectors" : 24,
    "aNoiseGain" : 0.1,
    "boundaryFalloffDist" : 0.015615,
    "nnInitialOffset" : 1.0,
    "ccInitialOffset" : 2.5,
    "overwrite_logs" : true,
    "skipMorph" : false,
    "lPerturb" : true,
    "iter" : "0",
    "lminRadius" : false,
    "lPolygon" : false,
    "off" : 1,
    "plotEvery" : 999,
    "checkEvery" : 999,
    "saveplots" : true,
    "vidframes" : true,
    "win_width" : 2050,
    "win_height" : 1025,
    "fov" : 40,
    "numpoints" : 41,
    "diffTol": 1e-12,
    "lengthScale": 25.0,
    "exponent": -1000.0
}
EOF

# Success/completion
exit 0
