#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
//#include "display.cpp"
#include <morph/Config.h>
#include <morph/display.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
#include <chrono>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Gdisplay;
using morph::HexGrid;
using morph::HdfData;
using morph::Tools;
using namespace morph;
using namespace std;
using namespace std::chrono;

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    string jsonfile = argv[1];
    //open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
    double dt = conf.getDouble("dt",0.0001);
    double Dn = conf.getDouble("Dn",5.0);
    double Dchi = conf.getDouble("DChi",5.0);
    double Dc = conf.getDouble("Dc",1.5);
    int scale = conf.getInt("scale",8);
    double xspan = conf.getDouble("xspan",5.0);
    int numsteps = conf.getInt("numsteps",100);
    int numAdjust = conf.getInt("numAdjust",1000000);
    int numprint = conf.getInt("numprint",95);
    string logpath = conf.getString("logpath","../logs");
    double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
    double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
    int numSectors = conf.getInt("numsectors",12);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    double nnInitialOffset = conf.getDouble("nnInitialOffset",1.0);
    double ccInitialOffset = conf.getDouble("ccInitialOffset", 2.5);
    bool overwrite_logs = conf.getBool("overwrite_logs",1);
    bool skipMorph  = conf.getBool("skipMorph",1);
    bool lminradius = conf.getBool("lminradius",0);
    unsigned int off = conf.getUInt("off",1);
    int NUMPOINTS = conf.getInt("NUMPOINTS",1);
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;
    bool lGraphics = true;
    /*
     * NOW create a log directory if necessary, and exit on any
     * failures.
     */
    if (morph::Tools::dirExists (logpath) == false) {
        morph::Tools::createDir (logpath);
        if (morph::Tools::dirExists (logpath) == false) {
            cerr << "Failed to create the logpath directory "
                 << logpath << " which does not exist."<< endl;
            return 1;
        }
	}
     else {
        // Directory DOES exist. See if it contains a previous run and
        // exit without overwriting to avoid confusion.
        if (overwrite_logs == false
            && (morph::Tools::fileExists (logpath + "/params.json") == true
                || morph::Tools::fileExists (logpath + "/guidance.h5") == true
                || morph::Tools::fileExists (logpath + "/positions.h5") == true)) {
            cerr << "Seems like a previous simulation was logged in " << logpath << ".\n"
                 << "Please clean it out manually, choose another directory or set\n"
                 << "overwrite_logs to true in your parameters config JSON file." << endl;
            return 1;
        }
    }
    ofstream gfile ( logpath + "/edges.out");
    ofstream jfile ( logpath + "/results.txt");
	ofstream degfile1 (logpath + "/degree1.data");
	ofstream degfile2 (logpath + "/degree2.data");
	ofstream degfile3 (logpath + "/degree3.data");


// initialise DRegion class setting scale
    DRegion M(scale,xspan,logpath,NUMPOINTS);
    cout << "before dissect_boundary " << endl;
    double hexWidth = M.Hgrid->hexen.begin()->d/2.0;
// include the analysis methods
    Analysis L;

    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
    unsigned int seed = time(NULL);
    // A rando2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);
// section for solving on the circle Tessllation
    cout << "just before creating circle Tessellation Solver S" << endl;
    //Readjust Dn for a single region
    FLT radius = 1.0;
    pair<FLT,FLT> centroid(0.0,0.0);
    M.setCreg();
    cout << "just after clearing region boundary" << endl;
    vector<std::pair<FLT, FLT>> centroids;
    centroids = M.dissectBoundary(); //dissect region boundary
    M.setRadialSegments(); //set radial segments and vertices for each region

// initialise the fields
    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
        //WARNING next 2 lines are just a temporary fix
   // Lcontinue=true;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        std::string ccstr = "c";
        char * ccst = new char[ccstr.length()+1];
        std::strcpy(ccst,ccstr.c_str());
        std::string nstr = "n";
        char * nst = new char[nstr.length()+1];
        std::strcpy(nst,nstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        ginput.read_contained_vals(ccst,M.cc);
        ginput.read_contained_vals(nst,M.nn);
    }
    else
    {
        unsigned int seed;
        //milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
        //seed = static_cast<unsigned int> (ms1.count());
        seed = off;
        morph::RandUniform<double> ruf1(seed);
        for (auto h : M.Hgrid->hexen) {
        // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
        // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
        // normal value. Close to boundary, noise is less.
            double choice = morph::Tools::randDouble();
		    if (choice > 0.5)
			{
                M.nn[h.vi] = - ruf1.get() * aNoiseGain +nnInitialOffset;
                M.cc[h.vi] = - ruf1.get() * aNoiseGain + ccInitialOffset;
			}
		    else
			{
                M.nn[h.vi] = - ruf1.get() * aNoiseGain +nnInitialOffset;
                M.cc[h.vi] = - ruf1.get() * aNoiseGain + ccInitialOffset;
            }
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                //double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                double bSig = 1.0;
                M.nn[h.vi] = (M.nn[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                M.cc[h.vi] = (M.cc[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
            } //end of if on boundary distance
        }//end of loop over region
        usleep(1000);
    } //end of else on Lcontinue
//work out the cutting template.
    vector<vector<hexGeometry::lineSegment>> regionLineSegments;
    regionLineSegments.resize(NUMPOINTS);
    for (int j=0;j<NUMPOINTS;j++)
    {
        regionLineSegments[j] = M.polygonSides(j, true);
         unsigned int sideSize = regionLineSegments[j].size();
         cout << "sideSize " << sideSize << " region " << j << endl;
         int count = 0;
         for (auto& h : M.Hgrid->hexen) {
             for (unsigned int i=0; i<sideSize; i++) {
                 if (M.hGeo->hexIntersectLineSegment(regionLineSegments[j][i],h)){
                     h.setFlags(HEX_IS_REGION_BOUNDARY);
                     //M.regionBound[j].push_back(h);
                     count++;
                  }
              }
         }
     }
    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    double rhoInit = 3.0;
    //double rhoInit = 5.0;
    array<float,3> cl_a = morph::Gdisplay::getJetColorF (0.1);
    array<float,3> cl_c = morph::Gdisplay::getJetColorF (0.28);
    array<float,3> cl_b = morph::Gdisplay::getJetColorF (0.58);
    array<float,3> cl_d = morph::Gdisplay::getJetColorF (0.80);
    array<float,3> offset = {{0, 0, 0}};
    int boundaryCount = 0;
     morph::Gdisplay ndisp(900, 900, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
     //plot the tessellation
     ndisp.resetDisplay (fix, eye, rot);
     for (int j=0;j<NUMPOINTS;j++) {
         // plot stuff here
         int boundaryCount = 0;
         int internalCount = 0;
         int totalCount = 0;
         for (auto h : M.Hgrid->hexen) {
             if (M.Creg[h.vi] == 1) {
                 ndisp.drawHex (h.position(), (h.d/2.0f), cl_c);
                 boundaryCount++;
             }
         	 if (M.Creg[h.vi] > 1) {
                 ndisp.drawHex (h.position(), (h.d * 2.0f), cl_d);
             }
         	 if (M.Cnbr[h.vi] == -1) {
                 ndisp.drawHex (h.position(), (h.d / 2.0f), cl_b);
             }
         }
         cout << " region " << j << " has boundary length " << boundaryCount << endl;
     }//end of loop on NUMPOINTS to plot boundaries
     usleep (1000000);
     cout << "before redrawDisplay 2 " << endl;
     ndisp.redrawDisplay();
     cout << "after redrawDisplay 2" << endl;
     usleep (1000000); // one hundred seconds
     ndisp.saveImage(logpath + "/Tessellation2.png");
     ndisp.closeDisplay();
    cerr << "d/2: " << hexWidth << endl;

// section for solving on the circle Tessllation
// if (skipMorph) return 0;
    cout << "just before creating circle Tessellation" <<endl;
    ksSolver S(scale, xspan, logpath, radius, centroid);
    FLT choice = ruf.get();
    if (choice > 0.5)
    {
        S.NN[h.vi] = - ruf1.get() * aNoiseGain +nnInitialOffset;
        S.CC[h.vi] = - ruf1.get() * aNoiseGain + ccInitialOffset;
    }
    else
    {
        S.NN[h.vi] =  ruf1.get() * aNoiseGain +nnInitialOffset;
        S.CC[h.vi] =  ruf1.get() * aNoiseGain + ccInitialOffset;
    }

    cout << "before time-stepping loop" << endl;
    for (int i=0;i<numsteps;i++) {
        S.step(dt, Dn, Dchi, Dc);
        if (i%numprint == 0 && lGraphics) {
            morph::Gdisplay disp(900, 900, 0, 0, "first run", rhoInit, 0.0, 0.0);
            disp.resetDisplay (fix, eye, rot);
            cout << "in print routine"<<endl;
            vector<double> normalnn;
	    normalnn.resize(M.n);
	    int countHex = 0;
            int regsize = M.n;

	    cout << "in loop over regions " << j << " size is " << regsize << endl;
            countHex += regsize;
	    vector<FLT> regionnn(regsize);
	    vector<FLT> tempnn(regsize);
	    vector<int> regionIdx(regsize);
	    int k = 0;
	    for (auto h : M.regionHex[j]){
                int index = h.vi;
                tempnn[k] = M.nn[index];
                regionIdx[k] = index;
		k++;
	    }
//cout << "normalise over the region then write  to normalised nn over hexGrid" << endl;
            regionnn = L.normalise(tempnn);
            for (unsigned int k=0;k<regsize;k++) {
                normalnn[regionIdx[k]] = regionnn[k];
            }
            std::array<float,3> pos = {M.centres[j].first, M.centres[j].second, 0};
            disp.drawHex (pos,offset,4.0*hexWidth,cl_d);
	    cout << "total number of hexes counted " << countHex << endl;
            for (auto h : M.Hgrid->hexen) {
            // cout << "just before drawHex  "<< h.vi << "normalnn " << normalnn[h.vi] << endl;
                if (M.Creg[h.vi] == 0) {
                    array<float,3> colour = morph::Gdisplay::getJetColorF(normalnn[h.vi]);
                    disp.drawHex(h.position(),hexWidth,colour);
                }
		        else if (M.Cnbr[h.vi] == -1 ) {
                    disp.drawHex(h.position(),hexWidth,cl_b);
                }
                else if (M.Creg[h.vi] == 1) {
                    disp.drawHex(h.position(), hexWidth, cl_c);
                }
                else if (M.Creg[h.vi] > 1) {
                    disp.drawHex(h.position(), 4*hexWidth, cl_d);
                }
            }
            cout << "just before redraw display 0" << endl;
            usleep (1000000);
            disp.redrawDisplay();
            usleep (1000000); // one hundred seconds
	     disp.saveImage(logpath + "/nnField.png");
	     disp.closeDisplay();
         } // end of print on numprint
     } //end of numsteps loop
//cout << " just after time step i = " << i << endl;

//code run at end of timestepping
//first save the  ofstream outFile;
         morph::HdfData data(fname);
         data.add_contained_vals("c",S.CC); //chemoattactant
         data.add_contained_vals("n",S.NN); //neuronal density
         data.add_contained_vals("x",M.Hgrid->d_x); //x coord of hex centres
         data.add_contained_vals("y",M.Hgrid->d_y); //y coord of hex centres
         data.add_val ("/Dchi", Dchi); //Chi parameter
         data.add_val ("/Dn", Dn); //Dn paramater
         data.add_val ("/Dc",Dc); //Dc parameter:
         cout << " just after writing data "  << endl;
// post run analysis
         vector <int> radiusDVector;
         vector <int> angleDVector;
	     vector <double> angleVector;
         vector <double> radiusVector;
         int degreeRadius;
         int degreeAngle;
	 int tempArea = 0;
	 int tempPerimeter = 0;
	 int angleOffset = 0;
	 int radiusOffset = 0;
	  double avDegreeAngle = 0;
	  double avDegreeRadius = 0;
	  double occupancy = 0;
          tempArea = 0;
	  tempPerimeter = 0;
	  int countRegions = 0;
      for (int j=0;j<NUMPOINTS;j++) {
	      if (M.regArea(j) != 0){
	          countRegions++;
              cout << "just before regNNFrac" << endl;
		      occupancy += M.regnnfrac(j);
              cout << "just after regNNFrac" << endl;
              tempArea = M.regArea(j);
              tempPerimeter = M.regPerimeter(j);

              angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, M.nnRegional(j));
              cout << "just ager sectorize Dangle" << endl;
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, M.nnRegional(j));
              cout << "just ager sectorize Dradius" << endl;
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              cout << "just before writing results first time" << endl;
              degfile1 << degreeAngle/2 << " " << degreeRadius << " " << occupancy << " " << tempArea<< " "<< tempPerimeter<< " xval " << M.centres[j].first << " yval " << M.centres[j].second <<endl<<flush;
              cout << "just after writing results first time" << endl;
		  }
	  } //end of loop on NUMPOINTs
      avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
      avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	  occupancy = occupancy / (1.0 * countRegions);
	  jfile <<Dn<<" "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy <<endl;
//end of integration on the original tesselation
    return 0;
};
