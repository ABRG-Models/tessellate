

#!/bin/zsh

# Set the parameters for displayCircles, to apply the cookie cutter method
# and compare it with the solutions that evolve in polygons

# Check we're running from BooleanDiffusion dir
CURDIR=$(pwd | awk -F '/' '{print $NF}')

if [ ! $CURDIR = "CookieCutter" ]; then
    echo "you neeed to run from sim/CookieCutter"
    exit 1
fi


            cat > ./dCircles.json <<EOF
{
    "logpath" : "./logsDisplayCircles",
    "dt": 0.0001,
    "Dn": $1,
    "Dchi" : $2,
    "Dc" : $3,
    "scale": 8,
    "xspan": 5.0,
    "numsteps": 5000,
    "numAdjust" : 100000,
    "numprint" : 4999,
    "Lcontinue" : 0,
    "LfixedSeed" : 0,
    "Lgraphics" : 1,
    "LDn" : 1,
    "numSectors" : 12,
    "aNoiseGain" : 0.1,
    "boundaryFalloffDist" : 0.0078,
    "nnInitialOffset" : 1.0,
    "overwrite_logs" : true,
    "skipMorph" : 0,
    "lPerturb" : 1,
    "lminRadius" : 0,
    "off" : 1,
    "NUMPOINTS" : 41,
    "lZero" : 1
}
EOF

# Success/completion
exit 0
