#include<morph/HexGrid.h>
#include<morph/ReadCurves.h>
#include<morph/HdfData.h>
#include <morph/tools.h>
#include <morph/ColourMap.h>
#include <utility>
#include <iostream>
#include <morph/ReadCurves.h>
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <random>
#include <algorithm>
#include <hdf5.h>
#include <unistd.h>
#include <bits/stdc++.h>
#include <sys/stat.h>
#include <sys/types.h>
using namespace std;


int main(int argc, char **argv) {
    int scale = std::stoi(argv[1]);
    FLT xspan = std::stof(argv[2]);
    std::string logpath = argv[3];
    FLT ds = 1.0/pow(2.0,scale-1);
    morph::HexGrid* Hgrid;
    FLT pspan = xspan/3.0; //make the parallelogram smaller than the original hexagon
    Hgrid = new morph::HexGrid(ds, xspan, 0.0, morph::HexDomainShape::Parallelogram);
    int rextent =  floor(0.5 * pspan / ds) + 1;
    int gextent = rextent;
    std::cout << " span of rextent " << rextent << std::endl;
    Hgrid->setParallelogramBoundary(rextent,gextent);
    std::cout << "after setParallelogramBoundary " << std::endl;
    std::cout << "length of rows " << Hgrid->d_rowlen << " number of rows " << Hgrid->d_numrows << std::endl;
    std::cout << " size of HexGrid " << Hgrid->d_size << std::endl;

    // Create a HexGrid Visual
    morph::Visual v(1600, 1000, "HexGrid");
    v.lightingEffects();
    morph::Vector<float, 3> offset = { 0.0f, -0.0f, 0.0f };
    morph::HexGridVisual<float>* hgv = new morph::HexGridVisual<float>(v.shaderprog, v.tshaderprog, Hgrid, offset);
    cout << "after morph::HexGridVisual" << endl;
    // Set up data for the HexGridVisual and colour hexes according to their state as being boundary/inside/domain, etc
    std::vector<float> colours (Hgrid->num(), 0.0f);
    static constexpr float cl_left = 1.0f;
    static constexpr float cl_right = 0.7f;
    static constexpr float cl_bottom = 0.3f;
    static constexpr float cl_top = 0.0f;
    static constexpr float cl_domain = 0.5f;
    // Note, HexGridVisual uses d_x and d_y vectors, so set colours according to d_flags vector
    for (unsigned int i = 0; i < Hgrid->num(); ++i) {
        if (i%301 == 0) {
            colours[i] = cl_left;
        }else if (i%301 == 300) {
            colours[i] = cl_right;
        }
        else if (i<301) {
            colours[i] = cl_bottom;
        }else if (i>= 58394 && i<58695) {
            colours[i] = cl_top;
        }else {
            colours[i] = cl_domain;
        }
    }
    hgv->cm.setType (morph::ColourMapType::Jet);
    hgv->zScale.setParams (0,0); // makes the output flat in z direction, but you still get the colours
    hgv->setScalarData (&colours);
    hgv->hexVisMode = morph::HexVisMode::HexInterp; // Or morph::HexVisMode::Triangles for a smoother surface plot
    hgv->finalize();
    v.addVisualModel (hgv);


    while (v.readyToFinish == false) {
        glfwWaitEventsTimeout (0.018);
        v.render();
    }

    return 0;
}

