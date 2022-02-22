

#!/bin/zsh

# Create json file for sHreal.cpp


            cat > ./sHreal.json <<EOF
{
    "dt" : 0.0000125,
    "epsilon" : 2.0,
    "g" : 1.02,
    "scale" : 8,
    "xspan" : 5.0,
    "numsteps" : 300000,
    "numAdjust" : 1000000,
    "numprint" : 19995,
    "logpath" : "./logsSwiftHohenberg",
    "boundaryFalloffDist" : 0.0078,
    "aNoiseGain" : 0.1,
    "numsectors" : 12,
    "Lcontinue" : false,
    "LfixedSeed" : false,
    "nnInitialOffset" : 1.0,
    "ccInitialOffset" : 2.5,
    "overwrite_logs" : 1,
    "skipMorph" : 1,
    "lminradius" : 0,
    "off" : 1,
    "plotevery" : 999,
    "saveplots" : true,
    "vidframes" : false,
    "win_width" : 1025,
    "NUMPOINTS" : 5
}
EOF

# Success/completion
exit 0
