#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <morph/Config.h>
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

/*
using morph::HexGrid;
using morph::HdfData;:
using morph::Tools;
using morph::Display
*/
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
    string iter = argv[2];
    //open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
    double dt = conf.getDouble("dt",0.0001);
    double Dn = conf.getDouble("Dn",0.0);
    double Dchi = conf.getDouble("Dchi",0.0);
    double Dc = conf.getDouble("Dc",0.0);
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
    bool Lgraphics = conf.getBool("Lgraphics", 1);
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    double nnInitialOffset = conf.getDouble("nnInitialOffset",1.0);
    double ccInitialOffset = conf.getDouble("ccInitialOffset", 2.5);
    bool overwrite_logs = conf.getBool("overwrite_logs",1);
    bool skipMorph  = conf.getBool("skipMorph",1);
    bool lminRadius = conf.getBool("lminRadius",0);
    unsigned int off = conf.getInt("off",1);
    int NUMPOINTS=conf.getInt("NUMPOINTS", 41);
    bool lZero = conf.getBool("lZero", 1);

    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;


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
    DRegion M(scale, xspan, logpath, NUMPOINTS);
    M.setCreg();
    M.setInternalBoundary();
    cout << "before dissect_boundary " << endl;
    vector<std::pair<FLT,FLT>> centroids;
    centroids = M.dissectBoundary(); //dissect region boundary
    M.setRadialSegments(); //set radial segments and vertices for each region
    int inReg = 0;
    inReg = M.setInnerRegion(); //set mask array for inner regions
    int rC = 0;
    //check for correct number of inner regions
    for (unsigned int j=0; j<NUMPOINTS;j++) {
        if (!M.innerRegion[j]) rC++;
    }
    if (inReg != rC) {
        cout << "Error: setInnerRegion returns " << inReg << " but outerRegions " << rC << endl;
        return -1;
    }
    else {
        cout << "Success: setInnerRegion no of outer regions " << inReg << endl;
    }
    cout << "after first setRadialSegments " << endl;
    FLT hexWidth = M.Hgrid->hexen.begin()->d/2.0;
// include the analysis methods
    Analysis L;

    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
    unsigned int seed = time(NULL);
    // A rando2yym uniform generator returning real/FLTing point types
    morph::RandUniform<FLT> ruf(seed);
// initialise with random field
    if (Lcontinue) {
	    morph::HdfData input (fname,1);
	    cout<< "just after first data read fname "<< fname << endl;
	    input.read_contained_vals("n",M.nn);
	    input.read_contained_vals("c",M.cc);
	    cout<< "just after input of nn and cc1"<< endl;
    }
    else {
		for (auto h : M.Hgrid->hexen) {
		    FLT choice = ruf.get();
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		    if (choice > 0.5)
			{
                M.nn[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                M.cc[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
			}
			else
			{
                M.nn[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                M.cc[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
			}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                M.nn[h.vi] = (M.nn[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                M.cc[h.vi] = (M.cc[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		    } //end of if on boundary distance
	    }//end of loop over region
    } //end of else on Lcontinue
    cout <<  "just after field creation" << endl;
    vector<double> fix(3.0, 0.0);
    vector<double> eye(3.0, 0.0);
    vector<double> rot(3.0, 0.0);
    double rhoInit = 3.0;
    //FLT rhoInit = 5.0;
    array<float,3> cl_a = morph::Gdisplay::getJetColorF (0.78);
    array<float,3> cl_c = morph::Gdisplay::getJetColorF (0.28);
    array<float,3> cl_b = morph::Gdisplay::getJetColorF (0.58);
    array<float,3> cl_d = {{0, 0, 0.0}};
    array<float,3> offset = {{0, 0, 0}};
    int boundaryCount = 0;
    unsigned int window = 2160;
    if (Lgraphics) {
        morph::Gdisplay isp(window, window, 0, 0, "A boundary", rhoInit, 0.0, 0.0);
     	cout << "after setting display" << endl;
     	isp.resetDisplay (fix, eye, rot);
    	cout << "after setting display" << endl;
    	//plot stuff here.
     	for (auto h : M.Hgrid->hexen) {
        	 if (M.Creg[h.vi] > 1) {
              	isp.drawHex (h.position(), (h.d/2.0f), cl_c);
         	}
         	else if (M.Creg[h.vi] == 1) {
                    isp.drawHex (h.position(), (h.d/2.0f), cl_d);
                 }
       //  else
        // {
          //   isp.drawHex (h.position(), offset, (h.d/2.0f), cl_b);
         //}
         }
     	for (int regcnt = 0; regcnt < NUMPOINTS;regcnt++) {
     	    std::array<FLT,3> pos = {M.centres[regcnt].first, M.centres[regcnt].second, 0};
    	    cout << " drawing centre of region x " << M.centres[regcnt].first << " y " << M.centres[regcnt].second << endl;
        	isp.drawHex (pos,offset,4.0*hexWidth,cl_d);
     	}
     	cout << "boundaryCount "<<boundaryCount<<endl;
     	usleep (1000000); // ten seconds
     	isp.redrawDisplay();
    	usleep (1000000); // ten seconds
     	isp.saveImage(logpath + "/Tessellation0" + iter + ".png");
     	isp.closeDisplay();
    }

    cerr << "d/2: " << hexWidth << endl;
    cout << "before time-stepping loop" << endl;
    for (int i=0;i<numsteps;i++) {
        M.step(dt, Dn, Dchi, Dc);
        if (i%numprint == 0 && Lgraphics) {
            morph::Gdisplay disp(window, window, 0, 0, "first run", rhoInit, 0.0, 0.0);
            disp.resetDisplay (fix, eye, rot);
            cout << "in print routine"<<endl;
            vector<FLT> normalnn;
	        normalnn.resize(M.n);
	        int countHex = 0;
	        for (int j=0;j<NUMPOINTS;j++) {
                if (M.innerRegion[j]) {
                unsigned int regsize = M.regionHex[j].size();
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
               // std::array<FLT,3> pos = {M.centres[j].first, M.centres[j].second, 0};
               // disp.drawHex (pos,offset,4.0*hexWidth,cl_d);
                }//end loop on inner regions
            } //end of loop over regions
	        cout << "total number of hexes counted " << countHex << endl;
            for (auto h : M.Hgrid->hexen) {
            // cout << "just before drawHex  "<< h.vi << "normalnn " << normalnn[h.vi] << endl;
                //if (M.Creg[h.vi] == 0) {
                array<FLT,3> colour = morph::Gdisplay::getJetColorF(normalnn[h.vi]);
	            disp.drawHex(h.position(),hexWidth,colour);
                //}
		        //else {
                //    disp.drawHex(h.position(),hexWidth,cl_d);
               // }
            }
            cout << "just before redraw display 0" << endl;
            usleep (1000000);
            disp.redrawDisplay();
            usleep (1000000); // one hundred seconds
	     disp.saveImage(logpath + "/nnField" + iter + ".png");
	     disp.closeDisplay();
         } // end of print on numprint
     } //end of numsteps loop
//cout << " just after time step i = " << i << endl;

//code run at end of timestepping
//first save the  ofstream outFile;
         morph::HdfData data(fname);
         data.add_contained_vals("c",M.cc); //chemoattactant
         data.add_contained_vals("n",M.nn); //neuronal density
         data.add_contained_vals("x",M.Hgrid->d_x); //x coord of hex centres
         data.add_contained_vals("y",M.Hgrid->d_y); //y coord of hex centres
         data.add_val ("/Dchi", Dchi); //Chi parameter
         data.add_val ("/Dn", Dn); //Dn paramater
         data.add_val ("/Dc",Dc); //Dc parameter:

         cout << " just after writing data "  << endl;
// post run analysis
// look at correlation between adjacent edges

        cout << "before correlate_edges" << endl;

        FLT avAbsCorrelation = M.correlate_edges(0,lZero);
        cout << "before correlate_edges zero morph" << endl;
        cout << " Average correlation whole region " << avAbsCorrelation << endl;
// look at correlation between random edges
         const int max_comp = NUMPOINTS*5;
         cout << "before random_correlate_edges" << endl;
         M.random_correlate(max_comp, 0, false, lZero);
         cout << "after random_correlate_edges" << endl;
//  cout<<"after correlate_edges" << endl;
         vector <int> radiusDVector;
         vector <int> angleDVector;
	     vector <FLT> angleVector;
         vector <FLT> radiusVector;
         int degreeRadius;
         int degreeAngle;
	     FLT tempArea = 0;
	     FLT tempPerimeter = 0;
		 int angleOffset = 0;
		 int radiusOffset = 0;
         for (int j=0;j<NUMPOINTS;j++) {
             vector<FLT> nnField = M.nnRegional(j);
             //if (M.innerRegion[j]){
			     int regionCount = 0;
                 gfile<<"in the degree loop" << endl;
                 tempArea = M.regArea(j); //area by hexes
                 tempPerimeter = M.regPerimeter(j); //perimeter by hexes

                  // digital version
                 angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, nnField);
                 degreeAngle = L.find_zeroDAngle(angleDVector);
                 gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                 angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, nnField);
                 angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M.find_max(angleVector,3);
                 degreeAngle = L.find_zeroAngle(angleVector,3);
                 gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                 radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, nnField);
                 degreeRadius = L.find_zeroDRadius(radiusDVector);
                 gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
                 degreeRadius = -100;
                 int newdegreeRadius = 0;
                 for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
			          radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, nnField);
			          newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			          if (newdegreeRadius > degreeRadius)
				          degreeRadius = newdegreeRadius;
		         }


                 gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


                 regionCount++;

	      //end of if on non-zero regions
      } //end of loop on NUMPOINTs
//writing out of the image file]
      cout << "after writing edges.out zero morph" << endl;

	  FLT avDegreeAngle = 0;
	  FLT avDegreeRadius = 0;
	  FLT occupancy = 0;
      tempArea = 0;
	  tempPerimeter = 0;
	  int countRegions = 0;
      for (int j=0;j<NUMPOINTS;j++) {
	      //if (M.innerRegion[j]){
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
              //end of loop on inner regions
	  } //end of loop on NUMPOINTs
      avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
      avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	  occupancy = occupancy / (1.0 * countRegions);
	  jfile <<Dn<<" "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<<endl;
//end of integration on the original tesselation



// section for solving on the circle Tessllation
// if (skipMorph) return 0;
    cout << "just before creating circle Tessellation" <<endl;
    vector<ksSolver> S;

    //set the vector for perturbation of the centres
    vector<FLT> rad;
    rad.resize(0);
    vector<FLT> maxRadius;
    vector<FLT> minRadius;
    maxRadius.resize(0);
    minRadius.resize(0);
    for (int j=0; j<NUMPOINTS; j++) {
        maxRadius.push_back(M.max_radius(j,true));
        minRadius.push_back(M.min_radius(j,true));
    }

    for (int j=0; j<NUMPOINTS; j++) {
        rad.push_back(minRadius[j]*3.0*ruf.get());
    }

    //set radius for creating circular regions also calculate the adjusted Dn values
    vector<FLT> DnVal;
    DnVal.resize(NUMPOINTS,0.0);
    vector<FLT> DchiVal;
    DchiVal.resize(NUMPOINTS,0.0);
    vector<FLT> DcVal;
    DcVal.resize(NUMPOINTS,0.0);
    vector<pair<FLT,FLT>> perturbCentres;
	for (int j = 0;j<NUMPOINTS;j++) {
        FLT angle = 2.0 * PI * ruf.get();
        std::pair<FLT, FLT> centroid =  M.baryCentre(j);
        centroid.first += rad[j]*cos(angle);
        centroid.second += rad[j]*sin(angle);
        perturbCentres.push_back(centroid);
    }
    for (int j=0; j<NUMPOINTS; j++) {
        FLT radius = 0;
        if (lminRadius) {
            radius = minRadius[j];
            DnVal[j] = Dn;
            DchiVal[j] = Dchi;
            DcVal[j] = Dc;
        }
        else {
            radius = maxRadius[j] + 3.0 * minRadius[j];
            FLT area = M.hexArea*M.regArea(j);
            DnVal[j] = Dn * sqrt(PI * radius * radius / area);
            DchiVal[j] = Dchi * sqrt(PI * radius * radius / area);
            DcVal[j] = Dc * sqrt(PI * radius * radius / area);
            cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
        }

        cout << " radius = " << radius << endl;
        S.push_back(ksSolver(scale, xspan, logpath, radius, perturbCentres[j]));
        cout << "in the loop populating the ksVector"<< j << " centroid.x " << perturbCentres[j].first - M.centres[j].first << " centroid.y " << perturbCentres[j].second - M.centres[j].second << endl;

    } //end of loop over NUMPOINTS creating the vector of hexGrids
	cout << "just after populating the ksVector"<<endl;
    // repopulate the regions and their boundaries with hex indices
    for (int j=0;j<NUMPOINTS;j++)
    {
        M.renewRegion(j,S[j].Hgrid->hexen);
        int rsize = M.regionHex[j].size();
        int ssize = M.regionBound[j].size();
        std::pair bcent = M.baryCentre(j);
        cout << "barycentre of " << j << " x " << bcent.first << " y " << bcent.second << " region size " << rsize << " bound size " << ssize << endl;
    }
    for (int j=0;j<NUMPOINTS;j++) {
        M.renewBoundary(j,S[j].Hgrid->hexen);
        int rsize = M.regionHex[j].size();
        int ssize = M.regionBound[j].size();
        std::pair bcent = M.baryCentre(j);
        cout << "barycentre of " << j << " x " << bcent.first << " y " << bcent.second << " region size " << rsize << " bound size " << ssize << endl;
     }

    // test the vertices are ordered correctly by increasing angle from the centroid
    for (int j=0; j<NUMPOINTS;j++) {
        M.testRegionVertices(j);
    }

// now draw the tessellation of overlapping circles
    if (Lgraphics) {
    	morph::Gdisplay mdisp(window, window, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
    	mdisp.resetDisplay (fix, eye, rot);
        hexWidth = M.Hgrid->hexen.begin()->d/2.0;
        for (int regcnt = 0; regcnt < NUMPOINTS;regcnt++) {
            std::array<FLT,3> pos = {M.centres[regcnt].first, M.centres[regcnt].second, 0};
            std::array<FLT,3> pos1 = {M.baryCentre(regcnt).first,M.baryCentre(regcnt).second,0};
            std::array<FLT,3> pos2 = {perturbCentres[regcnt].first,perturbCentres[regcnt].second,0};
            cout << " drawing centre of region x " << perturbCentres[regcnt].first << " y " << perturbCentres[regcnt].second << endl;
           // mdisp.drawHex (pos,offset,4.0*hexWidth,cl_a);
            mdisp.drawHex (pos1,offset,4.0*hexWidth,cl_b);
            mdisp.drawHex (pos2,offset,4.0*hexWidth,cl_d);
        }
        for (int j=0;j<NUMPOINTS;j++) {
            if (M.regArea(j) == 0) break;
// plot stuff here.
            int boundaryCount = 0;
            int internalCount = 0;
            cout << "size of Hexgrid j "<< j << " is " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.boundaryHex()) {
                    mdisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                    boundaryCount++;
                }
            }
        } //end of for loop on j
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        mdisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (100000); // one hundred seconds
        mdisp.saveImage(logpath + "/Tessellation1" + iter + ".png");
        mdisp.closeDisplay();
    }

// initialise the fields
    string gname = logpath + "/second.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(gname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        for (unsigned int j=0;j<NUMPOINTS;j++)
		{
            std::string ccstr = "c" + to_string(j);
		    cout << " j string " << to_string(j) << " length" << ccstr.length()<< endl;
			char * ccst = new char[ccstr.length()+1];
			std::strcpy(ccst,ccstr.c_str());
		    std::string nstr = "n" + to_string(j);
			char * nst = new char[nstr.length()+1];
			std::strcpy(nst,nstr.c_str());
			cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
	        ginput.read_contained_vals(ccst,S[j].CC);
	        ginput.read_contained_vals(nst,S[j].NN);
		}
    }
    else
    {
        for (unsigned int j=0;j<NUMPOINTS;j++)
		{
            unsigned int seed;
            milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
            seed = static_cast<unsigned int> (ms1.count());
            //seed = j + off;
            morph::RandUniform<FLT> ruf1(seed);
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                    FLT choice = ruf.get();
		    if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf1.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf1.get() * aNoiseGain + ccInitialOffset;
				}
		    else
				{
                    S[j].NN[h.vi] =  ruf1.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] =  ruf1.get() * aNoiseGain + ccInitialOffset;
				}
                /*
                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                    S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                    S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
                } //end of if on boundary distance
                */
            }//end of loop over region
            usleep(1000);
        }//end of loop over all regions
    } //end of else on Lcontinue
// begin second time stepping loop
    for (int i=0;i<numsteps;i++) {
        for (int j = 0;j<NUMPOINTS;j++) //loop over regions
        {
            S[j].step(dt, DnVal[j], DchiVal[j], DcVal[j]);
        }
    }
    //now populate regions with the final field values
    for (int j=0; j<NUMPOINTS; j++) {
        M.NN[j] = S[j].NN;
    }
    //sort over the circle boundary
    for (int j=0; j<NUMPOINTS; j++) {
        M.sortRegionBoundary(j);
    }
    //rotate each region by a random value from [0,2*PI]
    /*
    for (int j=0; j<NUMPOINTS; j++) {
        FLT phaseShift = 2.0 * PI * ruf.get();
        M.rotateRegion(j, phaseShift);
    }
    */
    //and find the phi value for zeros on the boundary
    //JMB 03/01/2020 must add code to print out the phi values of the zeros for each region
    //they are all appended in a file that can then be analysed
    string fileString = (logpath + "/zeroIndices.data");
    string fileString1 = (logpath + "/zeroPhi.data");
    for (int j=0;j<NUMPOINTS;j++) {
        vector<int> zeroIndices;
        vector<FLT> zeroPhi;
        vector<FLT> normalNN = M.meanzero_vector(M.sortedBoundaryNN[j]);
        //M.printFLTVect(fileString1,normalNN);
        cout << " sortedBundary size " << M.sortedBoundary[j].size() << " normalNN size " << normalNN.size() << endl;
        zeroIndices = L.find_zeroIndices(normalNN);
        cout << "zeroIndices size" << zeroIndices.size() << endl;
        int zsize = zeroIndices.size();
        for (int i=0; i<zsize; i++) {
            zeroPhi.push_back(M.sortedBoundaryPhi[j][zeroIndices[i]]);
        }
        M.printIntVect(fileString, zeroIndices);
        M.printFLTVect(fileString1, zeroPhi);
       // M.printIntVect(fileString,M.sortedBoundary[j]);
    }
    //set the cutter outline
    vector<vector<hexGeometry::lineSegment>> regionLineSegments;
    regionLineSegments.resize(NUMPOINTS);
    for (int j=0;j<NUMPOINTS;j++)
    {
       // M.clearRegionBound(j);
        regionLineSegments[j] = M.polygonSides(j, true);
        unsigned int sideSize = regionLineSegments[j].size();
        cout << "sideSize " << sideSize << " region " << j << endl;
        int count = 0;
        for (auto& h : S[j].Hgrid->hexen) {
            for (unsigned int i=0; i<sideSize; i++) {
                if (M.hGeo->hexIntersectLineSegment(regionLineSegments[j][i],h)){
                    h.setFlags(HEX_IS_REGION_BOUNDARY);
                    M.regionBound[j].push_back(h);
                    count++;
                }
            }
        }
        cout << "in region " << j << " size of region boundary " << M.regionBound[j].size() << " count " << count << endl;
    }

// we now draw the tessellation cut by the gingerbread cutter
     if (Lgraphics) {
        morph::Gdisplay ndisp(window, window, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
        ndisp.resetDisplay (fix, eye, rot);
        for (int j=0;j<NUMPOINTS;j++) {
            if (M.innerRegion[j]) {
            // plot stuff here
            int boundaryCount = 0;
            int internalCount = 0;
            int totalCount = 0;
            //std::array<FLT,3> pos1 = {M.baryCentre(j).first,M.baryCentre(j).second,0};
            std::array<FLT,3> pos2 = {perturbCentres[j].first,perturbCentres[j].second,0};
            cout << " drawing centre of region x " << perturbCentres[j].first << " y " << perturbCentres[j].second << endl;
            ndisp.drawHex (pos2,offset,4.0*hexWidth,cl_d);
            //cout << "size of region j "<< j << " is " << M.regionpHexList[j].size() << " size of Hexgrid " << S[j].Hgrid->num() << endl;
            //cout << "size of region j "<< j << " is " << S[j].Hgrid->num() << " size of Hexgrid " << S[j].Hgrid->num() << endl;
            for (auto h : S[j].Hgrid->hexen) {
                if (h.testFlags(HEX_IS_REGION_BOUNDARY) == true) {
                    ndisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                    boundaryCount++;
                }
            }
            cout << " region " << j << " has boundary length " << boundaryCount << endl;
            /*
            // Draw small hex at boundary centroid
            array<FLT,3> c;
            c[2] = 0;
            c[0] = S[j].Hgrid->boundaryCentroid.first;
            c[1] = S[j].Hgrid->boundaryCentroid.second;
            cout << "d/2: " << S[j].Hgrid->hexen.begin()->d/4.0f << endl;
            ndisp.drawHex (c,  hexWidth*2.0f, cl_c);
            cout << "boundaryCentroid x,y: " << c[0] << "," << c[1] << endl;
            */
            } //end of loop on inner regions
        }//end of loop on NUMPOINTS to plot boundaries
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        ndisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (1000000); // one hundred seconds
        ndisp.saveImage(logpath + "/Tessellation2" + iter + ".png");
        ndisp.closeDisplay();
    }
    // print out the final field as with the gingerbread cutter tessellation
    if (Lgraphics) {
        //set up display
            morph::Gdisplay mmdisp(window, window, 0, 0, "Boundary 2", rhoInit, 0.0, 0.0);
            mmdisp.resetDisplay (fix, eye, rot);
            for (int j = 0;j<NUMPOINTS;j++) {
               if (M.innerRegion[j]) {
                cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
                cout << " region " << j << " size is " << regsize << endl;
                vector<FLT> tempNN;
                vector<FLT> regionNN;
        	    for (auto h : S[j].Hgrid->hexen) {
                    if (M.hexInRegion(j,h)) {
                        tempNN.push_back(S[j].NN[h.vi]);
                        cout << "tempNN " << S[j].NN[h.vi] << " h.vi " << h.vi  << endl;
                    }
        	    }
        //normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(tempNN);
                int idx = 0;
        	    for (auto h : S[j].Hgrid->hexen)
        	    {
                     if (M.hexInRegion(j,h)) {
                //std::pair<FLT, FLT> inPoint;
                //inPoint.first = h.x;
                //inPoint.second = h.y;
                //cout << " inRegion " << " xval " << inPoint.first << " yval " << inPoint.second <<endl;
   			        array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionNN[idx]);
                    mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
                    cout << "just after drawHex 2  value " << regionNN[h.vi] << endl;
                        /*
                        std::array<FLT,3> pos2 = {perturbCentres[j].first,perturbCentres[j].second,0};
                        mmdisp.drawHex (pos2,offset,4.0*hexWidth,cl_d);
                        */
                        idx++;
                    }
                }
                cout << "idx in new Hex routine  " << idx << endl;
                /*
                for (auto h : S[j].Hgrid->hexen) {
                   // if (h.boundaryHex()) {
                   if (h.testFlags(HEX_IS_REGION_BOUNDARY) == true)
                        //array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionNN[h.vi]);
                        //mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
                        mmdisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                        boundaryCount++;
                    }
                */
                } //end of loop on inner regions
            } //end of loop over regions
            cout << "just before redraw display 1" << endl;
            mmdisp.redrawDisplay();
            cout << "just after redraw display 1" << endl;
            cout << "just after to_string"<<endl;
            //int step = i / numprint;
            //string str = std::to_string(step);
            mmdisp.saveImage(logpath + "/nnField2" + iter + ".png");
            usleep (1000000); // one hundred seconds
            cout << "just after saveImage 1" << endl;
            mmdisp.closeDisplay();
     } //end of printing of gingerbread cut field
    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData gdata(gname);
    for (unsigned int j=0;j<NUMPOINTS;j++) {
	std::string nstr = "n" + to_string(j);
	char * nst = new char[nstr.length()+1];
	//std::copy(nstr.begin(),nstr.end(),nst);
	std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
        char * ccst = new char[ccstr.length()+1];
        std::strcpy(ccst,ccstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        gdata.add_contained_vals(ccst,S[j].CC);
        gdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M.X[0]);
        //data.add_contained_vals("Y",M.X[1]);
     }
     gdata.add_val ("/Dchi", Dchi);
     gdata.add_val ("/Dn", Dn);
     gdata.add_val ("/Dc",Dc);

     cout << " just after writing data "  << endl;

      gfile << endl << "analysis on first morphing iteration " << endl;
    for (int j=0;j<NUMPOINTS;j++) {
        //if (M.innerRegion[j]){
            int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);
            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors,S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors,S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            // radiusVector = M.meanzero_vector(radiusVector);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;
            ///radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
                radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                newdegreeRadius = L.find_zeroRadius(radiusVector,3);
                if (newdegreeRadius > degreeRadius)
                    degreeRadius = newdegreeRadius;
                }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
            degfile1 << degreeAngle/2 << " " << degreeRadius << " " << occupancy << " " << tempArea<< " "<< tempPerimeter<< " xval " << M.centres[j].first << " yval " << M.centres[j].second <<endl<<flush;
            cout << "just after writing results first time" << endl;
         //end of if on non-zero regions
    } //end of loop on NUMPOINT

    if (skipMorph) return 0;
    // clear the edges
    M.edges_clear();
    // redissect the boundaries
    cout << "just before renewDissect morph2 " << endl;
    for (int j=0;j<NUMPOINTS;j++)
	{
        M.renewDissect(j, 1);
    }

    cout << "just after renewDissect morph2 " << endl;
//computing average values over all regions
    avAbsCorrelation = 0;
    tempPerimeter = 0;
    tempArea = 0;
    countRegions = 0;
    M.random_correlate(max_comp, 1, true, lZero);
    cout << "just after random correlate_edges morph2 " << endl;
    avAbsCorrelation = 0;
    cout << "just after M.correlate_edges, third pass" << endl;
    for (int j=0;j<NUMPOINTS;j++)
	{
//        if (M.innerRegion[j]) {
	    countRegions++;
        avAbsCorrelation += M.renewcorrelate_edges(j,1,lZero);
//        }
	} //end of loop on NUMPOINTs
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	// jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	cout << Dn <<"  "<< Dchi <<" "<< Dc <<" " << avAbsCorrelation << " after morphing  2" << endl;

        // write the edge vectors all interpolated to a uniform size
        std::map<int, vector<FLT>>::iterator ptr;
        std::string sideNN = logpath + "/edgeNN.data";
        int ecount = 0;
        for (ptr = M.edgeNN.begin(); ptr != M.edgeNN.end(); ptr++) {
            vector<FLT> tempVect;
            tempVect = M.normaliseVect(ptr->second);
            M.printFLTVect(sideNN, tempVect);
            cout << "edgeNN key " << ptr->first << " element " << ecount << " vector size " << ptr->second.size() << " edges size " << M.edges[ptr->first].size() <<  endl;
            ecount++;
        }
//end of integration after solving on polygonal regions via ksSolver
    return 0;
}
