#include <morph/tools.h>
#include <morph/HexGrid.h>
#include <morph/HdfData.h>
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
#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;

//using morph::HexGrid;
//using morph::HdfData;
//using morph::Tools;
//using morph::Display
using namespace morph;
using namespace std;

int main (int argc, char **argv)
{
    ksSolver S;

    pair<float,float> v1 = make_pair (-0.7f, -0.3f);
    pair<float,float> v2 = make_pair (0.4f, -0.2f);
    pair<float,float> v3 = make_pair (1.0f, 0.0f);
    pair<float,float> v4 = make_pair (-0.4f, 0.9165f);
    pair<float,float> v5 = make_pair (-0.3f, 0.1f);
    cout << "ater making pairs" << endl;;
    morph::BezCurve<float> c1(v1,v2);
    morph::BezCurve<float> c2(v2,v3);
    morph::BezCurve<float> c3(v3,v4);
    morph::BezCurve<float> c4(v4,v5);
    morph::BezCurve<float> c5(v5,v1);
    cout << "after making BezCurves" << endl;
    morph::BezCurvePath<float> bound;
    cout << "after making BezCurvePath" << endl;
    bound.addCurve(c1);
    cout << "after adding curve" << endl;
    bound.addCurve(c2);
    cout << "after adding curve" << endl;
    bound.addCurve(c3);
    cout << "after adding curve" << endl;
    bound.addCurve(c4);
    bound.addCurve(c5);
    cout << "after adding curve" << endl;
    cout << "after adding curve" << endl;
    float radius = 1.0;
    double xspan = 5.0;
    int scale = 8;
    std::string logpath = "./logs";
    std::pair<double,double> centroid (0.0, -0.0);
    S = ksSolver(scale, xspan, logpath, radius, centroid);// now draw the intial tesselation
    float rhoInit = 3.0;
    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    array<float,3> cl_a = morph::Tools::getJetColorF (0.78);
    array<float,3> cl_c = morph::Tools::getJetColorF (0.28);
    array<float,3> cl_b = morph::Tools::getJetColorF (0.58);
    array<float,3> cl_d = morph::Tools::getJetColorF (0.00);
    array<float,3> offset = {{0, 0, 0}};
    morph::Gdisplay mdisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
    mdisp.resetDisplay (fix, eye, rot);
    double hexWidth = S.Hgrid->hexen.begin()->d/2.0;
// plot stuff here.
    int internalCount = 0;
    int boundaryCount = 0;
    cout << "size of Hexgrid  "<< " is " << S.Hgrid->num() << endl;
    for (auto h : S.Hgrid->hexen) {
            if (h.boundaryHex()) {
                mdisp.drawHex (h.position(), (h.d/2.0f), cl_a);
    			boundaryCount++;
            }
        }
    usleep (1000000);
    cout << "before redrawDisplay 2 " << endl;
    mdisp.redrawDisplay();
    cout << "after redrawDisplay 2" << endl;
    usleep (100000); // one hundred seconds
    mdisp.saveImage(logpath + "/Tesselation0.png");
    mdisp.closeDisplay();
    morph::Gdisplay ndisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
    ndisp.resetDisplay (fix, eye, rot);
    boundaryCount = 0;
    internalCount = 0;
    int totalCount = 0;
    std::pair<float,float> regionCentroid;
    std::vector<std::list<morph::Hex>::iterator> regionpHexList; //vector of vectors of pHexes (iterators)
    regionpHexList = S.Hgrid->getRegion(bound, regionCentroid, false);
    cout << "size of region "<< " is " << regionpHexList.size() << " size of Hexgrid " << S.Hgrid->num() << endl;
    for (auto h : regionpHexList) {
        totalCount++;
        if (h->testFlags(HEX_IS_REGION_BOUNDARY) == true)
        {
             ndisp.drawHex (h->position(), (h->d/2.0f), cl_a);
             boundaryCount++;
        }
        else
        {
            ndisp.drawHex (h->position(), offset, (h->d/2.0f), cl_b);
            internalCount++;
        }
    }
    usleep (1000000);
    cout << "before redrawDisplay 2 " << endl;
    ndisp.redrawDisplay();
    cout << "after redrawDisplay 2" << endl;
    usleep (10000000); // one hundred seconds
    cout << "after redrawDisplay 2" << endl;
    ndisp.redrawDisplay();
    usleep (10000000); // one hundred seconds
    ndisp.saveImage(logpath + "/Tesselation1.png");
    cout << "after redrawDisplay 2" << endl;
    ndisp.closeDisplay();
        cout << "boundaryCount 2 "<<boundaryCount<< " internalCount 2 "<< internalCount << " totalCount 2 = " << totalCount << "  region size " << regionpHexList.size() << endl;

   return 0;
}
