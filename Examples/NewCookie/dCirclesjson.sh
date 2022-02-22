#!/bin/zsh

# Set the parameters for displayCircles, to apply the cookie cutter method
# and compare it with the solutions that evolve in polygons

# Check we're running from BooleanDiffusion dir
CURDIR=$(pwd | awk -F '/' '{print $NF}')

if [ ! $CURDIR = "NewCookie" ]; then
    echo "you neeed to run from NewCookie"
    exit 1
fi


            cat > ./dCircles.json <<EOF
{
    "logpath" : "./logsdCircles",
    "dt": 0.0001,
    "Dn": $1,
    "Dchi" : $2,
    "Dc" : $3,
    "scale": 9,
    "xspan": 5.0,
    "numsteps": 1000000,
    "numAdjust" : 100000,
    "numprint" : 1000,
    "Lcontinue" : false,
    "LfixedSeed" : false,
    "Lgraphics" : true,
    "LDn" : true,
    "numSectors" : 12,
    "aNoiseGain" : 1.0,
    "boundaryFalloffDist" : 0.0078,
    "nnInitialOffset" : 1.0,
    "ccInitialOffset" : 2.5,
    "overwrite_logs" : true,
    "skipMorph" : false,
    "lPerturb" : false,
    "LDn" : true,
    "off" : 1,
    "numpoints" : 41,
    "lZero" : true,
    "plotEvery" : 1000,
    "diffTol" : 1e-8
}
EOF

# Success/completion
exit 0
