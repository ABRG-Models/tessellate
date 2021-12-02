

#!/bin/zsh

# Explore the 2 gene model, initialised with a Gaussian hump in Gene a.


            cat > ./singleCircle.json <<EOF
{
    "dt" : 0.0001,
    "Dn" : 10.0,
    "DChi" : 5.0,
    "Dc" : 3.3,
    "scale" : 8,
    "xspan" : 5.0,
    "numsteps" : 1000,
    "numAdjust" : 1000000,
    "numprint" : 95,
    "logpath" : "../logs",
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
    "NUMPOINTS" : 5
}
EOF

# Success/completion
exit 0
