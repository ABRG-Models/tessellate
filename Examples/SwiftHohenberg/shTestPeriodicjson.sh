#!/bin/zsh

# Create json file for sHreal.cpp

            cat > ./shTestPeriodic.json <<EOF
{
    "dt" : 0.00005,
    "epsilon" : 0.1,
    "g" : 0.98,
    "k0" : 1.0,
    "scale" : 7,
    "xspan" : 6.0,
    "numsteps" : 100,
    "numAdjust" : 1000000,
    "numprint" : 99,
    "logpath" : "./logsSwiftHohenberg",
    "boundaryFalloffDist" : 0.0078,
    "aNoiseGain" : 0.1,
    "numsectors" : 12,
    "red" : 100,
    "green" : 160,
    "fov" : 70,
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
    "nonLocal" : 45,
    "saveplots" : true,
    "vidframes" : true,
    "showfft" : false,
    "win_width" : 2050,
    "NUMPOINTS" : 5,
    "nbins" : 100
}
EOF

# Success/completion
exit 0
