#!/bin/zsh

# Create json file for sHreal.cpp

            cat > ./shComplex.json <<EOF
{
    "dt" : 0.00001,
    "epsilon" : 0.1,
    "g" : 0.98,
    "k0" : 1.0,
    "scale" : 8,
    "xspan" : 6.0,
    "numsteps" : 100000,
    "numCheck" : 25000,
    "numprint" : 999,
    "logpath" : "./logsSwiftHohenberg",
    "boundaryFalloffDist" : 0.0078,
    "aNoiseGain" : 0.1,
    "numsectors" : 12,
    "red" : 100,
    "green" : 160,
    "fov" : 70,
    "Lcontinue" : true,
    "LfixedSeed" : false,
    "nnInitialOffset" : 1.0,
    "ccInitialOffset" : 2.5,
    "x_default": 0.0,
    "y_default": 0.0,
    "wratio": 0.867,
    "radius": 2.0,
    "ROIwid": 1.4,
    "lengthScale": 29.0,
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
    "nbins" : 100,
    "gaussBlur" : 1
}
EOF

# Success/completion
exit 0
