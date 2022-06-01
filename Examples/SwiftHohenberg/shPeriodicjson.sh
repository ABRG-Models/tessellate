#!/bin/zsh

# Create json file for sHreal.cpp

            cat > ./shPeriodic.json <<EOF
{
    "dt" : 0.0000125,
    "epsilon" : 2.0,
    "g" : 1.02,
    "scale" : 8,
    "xspan" : 3.0,
    "numsteps" : 4000000,
    "numAdjust" : 1000000,
    "numprint" :995,
    "logpath" : "./logsSwiftHohenberg",
    "boundaryFalloffDist" : 0.0078,
    "aNoiseGain" : 0.1,
    "numsectors" : 12,
    "red" : 100,
    "green" : 160,
    "fov" : 30,
    "Lcontinue" : false,
    "LfixedSeed" : false,
    "nnInitialOffset" : 1.0,
    "ccInitialOffset" : 2.5,
    "x_default": 0.0,
    "y_default": 0.0,
    "wratio": 0.867,
    "overwrite_logs" : 1,
    "skipMorph" : 1,
    "lminradius" : 0,
    "off" : 1,
    "plotevery" : 999,
    "saveplots" : true,
    "vidframes" : true,
    "win_width" : 2050,
    "NUMPOINTS" : 5
}
EOF

# Success/completion
exit 0
