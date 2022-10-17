#include <morph/tools.h>
#include <morph/HdfData.h>
#include <morph/Random.h>
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
#include <string>
#include <cctype>
#define PI 3.1415926535897932
#define NUMPOINTS 41 //just the A-E rows.
//#define NUMPOINTS 79 //just the A-E rows.

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HdfData;
using morph::Tools;
using namespace std;



int main (int argc, char **argv)
{
    FLT sqSide = stod(argv[1]);
    unsigned int nTicks = stoi(argv[2]);
    FLT stretch = stod(argv[3]);
    vector <std::pair <float, float>> centres; //seed points for chess squares
    unsigned int numpoints = nTicks*nTicks;
    centres.resize(numpoints);
    ofstream afile("./rectCentres.inp");
    ofstream bfile ("./rectInfo.txt");
    ofstream cfile ("./bezRect.h");
    int jump = (nTicks+1)/2;
    FLT overN = 1.0/nTicks;
    FLT overStretch = 1.0/stretch;
    FLT sqSidex = sqSide*stretch;
    FLT sqSidey = sqSide*overStretch;
    int count = 0;
    for (unsigned int i=0; i<nTicks; i++) {
        FLT x = 2.0*sqSide*(1.0*(i+1) - 1.0*jump)*overN*stretch;
        for (unsigned int j=0; j<nTicks; j++) {
            FLT y = 2.0*sqSide*(1.0*(j+1) -1.0*jump)*overN*overStretch;
            centres[count].first = x;
            centres[count].second = y;
            afile <<  centres[count].first << "  " << centres[count].second  << std::endl;
            count++;
        }
    }
    bfile << "sqSide = " << sqSide << " nTicks " << nTicks << " stretch " << stretch << " chessboard size " << count << std::endl;
    bfile << "sqSidex " << sqSidex << " sqSidey " << sqSidey << std::endl;

    char sideX[10];
    sprintf(sideX, "%f", sqSidex);
    char sideY[10];
    sprintf(sideY, "%f", sqSidey);
    cfile << "pair<float,float> v1 = make_pair (-" << sideX <<", -" << sideY << ");" << std::endl;
    cfile << "pair<float,float> v2 = make_pair (" << sideX << ", -" << sideY << ");" << std::endl;
    cfile << "pair<float,float> v3 = make_pair (" << sideX << ", " << sideY << ");" << std::endl;
    cfile << "pair<float,float> v4 = make_pair (-" << sideX << ", " << sideY << ");" << std::endl;
    cfile << "morph::BezCurve<float> c1(v1,v2);" << std::endl;
    cfile << "morph::BezCurve<float> c2(v2,v3);" << std::endl;
    cfile << "morph::BezCurve<float> c3(v3,v4);" << std::endl;
    cfile << "morph::BezCurve<float> c4(v4,v1);" << std::endl;
    cfile << "morph::BezCurvePath<float> bound;" << std::endl;
    cfile << "bound.addCurve(c1);" << std::endl;
    cfile << "bound.addCurve(c2);" << std::endl;
    cfile << "bound.addCurve(c3);" << std::endl;
    cfile << "bound.addCurve(c4);" << std::endl;

    return 0;
};
