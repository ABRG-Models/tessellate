

#!/bin/zsh

# Explore the 2 gene model, initialised with a Gaussian hump in Gene a.

# Check we're running from eanDiffusion dir
CURDIR=$(pwd | awk -F '/' '{pr $NF}')

if [ ! $CURDIR = "sim" ]; then
    echo "you neeed to run from sim"
    exit 1
fi


            cat > ./kst.json <<EOF
{
    "dt" : 0.0001,
    "Dn" : 12.0,
    "DChi" : 0.0,
    "Dc" : 3.6,
    "scale" : 8,
    "xspan" : 5.0,
    "numsteps" : 10000,
    "numAdjust" : 100000,
    "numprint" : 9995,
    "logpath" : "./logsKSTriangle",
    "Lgraphics" : true,
    "LDn" : true,
    "boundaryFalloffDist" : 0.0078,
    "aNoiseGain" : 0.1,
    "numsectors" : 12,
    "Lcontinue" : false,
    "nnInitialOffset" : 1.0,
    "ccInitialOffset" : 2.5,
    "overwrite_logs" : true,
    "skipMorph" : true,
    "lPerturb" : false,
    "LfixedSeed" : false,
    "ratio" : 1.0,
    "pRatio" : 0.0,
    "saveplots":true,
    "vidframes":true,
    "win_width" : 1025
}
EOF

# Success/completion
exit 0
