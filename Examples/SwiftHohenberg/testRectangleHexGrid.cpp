#include <morph/tools.h>
#include <morph/ReadCurves.h>
#include <morph/HdfData.h>
#include <morph/BezCurve.h>
#include <morph/BezCurvePath.h>
#include <morph/BezCoord.h>
#include <morph/HexGrid.h>
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
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#ifndef PI
#define PI 3.1415926535897932f
#endif

#ifndef NO_N
#define NO_N
#define nE(hi) (this->Hgrid->d_ne[hi])
#define HAS_nE(hi) (this->Hgrid->d_ne[hi] == -1 ? false : true)

#define nW(hi) (this->Hgrid->d_nw[hi])
#define HAS_nW(hi) (this->Hgrid->d_nw[hi] == -1 ? false : true)

#define nNE(hi) (this->Hgrid->d_nne[hi])
#define HAS_nNE(hi) (this->Hgrid->d_nne[hi] == -1 ? false : true)

#define nNW(hi) (this->Hgrid->d_nnw[hi])
#define HAS_nNW(hi) (this->Hgrid->d_nnw[hi] == -1 ? false : true)

#define nSE(hi) (this->Hgrid->d_nse[hi])
#define HAS_nSE(hi) (this->Hgrid->d_nse[hi] == -1 ? false : true)

#define nSW(hi) (this->Hgrid->d_nsw[hi])
#define HAS_nSW(hi) (this->Hgrid->d_nsw[hi] == -1 ? false : true)
#endif

int main(int argc, char **argv) {
    int scale = 4; //sets hext to hex distance
    FLT xspan = 6.0; //width of the hexGrid
    int n; //number of hexes in the HexGrid
    FLT ds; //distance used in numerical approximations
    FLT s = pow(2.0, scale-1);
    ds = 1.0/s;
    morph::HexGrid* Hgrid;
    FLT xwidth = 2.125; //width of the rectangle
    FLT ywidth = 2.125; //height of the rectangle
    Hgrid = new morph::HexGrid(ds, xspan, 0.0, morph::HexDomainShape::Rectangle);
    std::cout << "after setting rectangle boundary " << std::endl;
    Hgrid->setRectangularBoundary(xwidth, ywidth);
    std::cout << "after setting rectangle boundary " << std::endl;
    Hgrid->populate_d_neighbours();
    n = Hgrid->num();
    std::cout << "after populating d_neighbours n " << n << std::endl;
    std::cout << "after creating HexGrid ds =  " << ds << std::endl;
    std::cout << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << std::endl;
    std::cout << " rowlen " << Hgrid->d_rowlen << " numrows " << Hgrid->d_numrows << std::endl;
    int rowlen = Hgrid->d_rowlen;
    int numrows = Hgrid->d_numrows;
    FLT hexWidth = Hgrid->width();
    FLT hexDepth = Hgrid->depth();
    std::cout << "width " << hexWidth << std::endl <<  " depth " << hexDepth << std::endl;
    if (rowlen*numrows != n) {
        std::cerr << "error rowlen*numrows " << rowlen*numrows << " S.n " << n << std::endl;
        //return -1;
    }

    std::list<morph::Hex> boundaryHexes = Hgrid->getBoundary();
    for (auto h : boundaryHexes) {
        std::cout << " hex " << h.vi << " red " << h.ri << " green " << h.gi << " x " << h.x << " y " << h.y << std::endl;
    }

    return 0;
}
