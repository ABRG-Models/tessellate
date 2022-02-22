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
    double D_A = conf.getDouble("D_A",1.0);
    double D_B = conf.getDouble("D_B",20.0);
    double k1 = conf.getDouble("k1",0.01);
    double k2 = conf.getDouble("k2",1.0);
    double k3 = conf.getDouble("k3",1.0);
    double k4 = conf.getDouble("k4",1.7);
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
    bool LDn = conf.getBool("LDn", true);
    unsigned int off = conf.getInt("off",1);
    int NUMPOINTS=conf.getInt("NUMPOINTS", 41);
    bool lZero = conf.getBool("lZero", 1);
    double diffTol = conf.getDouble("diffTol", 0.000001);
    float lengthScale = conf.getFloat("lengthScale", 29.0);
    int checkEvery = conf.getInt("checkEvery", 1000);
    bool Lperturb = conf.getInt("Lperturb", true);

    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;
    std::cout << "D_A " << D_A << " D_B " << D_B << " k1 " << k1 << " diffTol " << diffTol << " NUMPOINTS " << NUMPOINTS << " Lperturb " << Lperturb << " LDn " << LDn << std::endl;
    std::cerr << "D_A " << D_A << " D_B " << D_B << " k1 " << k1 << " diffTol " << diffTol << " NUMPOINTS " << NUMPOINTS << " Lperturb " << Lperturb << " LDn " << LDn << std::endl;


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
   //M.setInternalBoundary(); do we still need this as well as setCreg
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
    milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    seed = static_cast<unsigned int> (ms1.count());
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
#ifdef COMPILE_PLOTTING
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
#endif

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate schSolver
     */
    vector<schSolver> S;
    /*
     * Set boundaries based on polygons derived from the Dregion tesselation.
     * stored in curvedBoundary
     */
    M.populateBoundPolygon(1);
    cout << "just after setting polygonal boundaries " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        S.push_back(schSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j],lengthScale));
        cout << "in the loop populating the schVector morph0 "<< j <<endl;
    }
    //now set the parameters for each solver
    for (int j=0; j<numpoints;j++) {
        S[j].setParams(D_A, D_B, k1, k2, k3, k4);
    }
    cout << "first setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
    // repopulate the regions
    cout << "just before populating the regions Morph 0" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        M.renewRegion(j,S[j].Hgrid->hexen);
    }
    cout << "just before populating the inner boundary Morph 0" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        M.renewBoundary(j,S[j].Hgrid->hexen);
        M.renewCentroids(j);
    }
    cout << "just before calculating regionSize " << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << " Region size " << M.regionHex[j].size() << endl;
    }
	cout << "just after populating the regions from the schSolver vector"<<endl;
    //clear global edges map
    M.edges_clear();
    cout << "just before renewDissect first time" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
        M.renewDissect(j,0);
    }
    cout << "Edges size " << M.edges.size() << endl;
/*
    FLT nndiff = 0.0;
    std::vector<FLT> nnbef;
    std::vector<FLT> nnnow;
    nnbef.resize(M.n,10000.0);
    nnnow.resize(M.n,0.0);

    //start of time-stepping loo
    cerr << "d/2: " << hexWidth << endl;
    cout << "before time-stepping loop" << endl;
    for (int i=0;i<numsteps;i++) {
        M.step(dt, Dn, Dchi, Dc);
*/
    // begin time stepping loop unmorphed grid solved via schSolver
    // set up vectors for determining convergence
    std::vector<FLT> Adiff;
    std::vector<std::vector<FLT>> Apre;
    std::vector<std::vector<FLT>> Acurr;
    Adiff.resize(numpoints);
    Apre.resize(numpoints);
    Acurr.resize(numpoints);
    FLT AdiffSum = 0.0f;
  //initilise all Apre vectors to above possible field
    if (!Lcontinue){
        for (int j=0; j<numpoints;j++) {
            Apre[j].resize(S[j].A.size(),10000.0);
        }
    }
    else {
        for (int j=0; j<numpoints; j++) {
            Apre[j] = S[j].A;
        }
    }

    // begin morph0 time stepping loop
    for (int i=0;i<numsteps;i++) {
   	for (unsigned int j = 0;j<numpoints;j++) { //loop over all regions, only step internal ones
            S[j].step(dt);
            cout <<"after step"<<endl;
            if (i%plotevery == 0) {
                AdiffSum = 0.0;
                Acurr[j] = S[j].A;
                Adiff[j] = L.normedDiff(Apre[j], Acurr[j]);
                Apre[j] = Acurr[j];
                AdiffSum += fabs(Adiff[j]);
                cout <<"after setting Apre etc.." << endl;
            }
        } //end of loop over regions
        if (i%plotevery == 0) {
            cerr << "AdiffSum " << AdiffSum << " i = " << i << endl;
            if (AdiffSum/(numpoints*1.0) < diffTol) {
                break;
            }
        }

#ifdef COMPILE_PLOTTING
        if (i%numprint == 0) {
            morph::Gdisplay disp(window, window, 0, 0, "first run", rhoInit, 0.0, 0.0);
            disp.resetDisplay (fix, eye, rot);
            cout << "in print routine"<<endl;
            int countHex = 0;
            for (int j=0;j<NUMPOINTS;j++) {
                if (M.innerRegion[j]) {
                    vector<FLT> regionA;
//cout << "normalise over the region then write  to normalised nn over hexGrid" << endl;
                    regionA = L.normalise(S[j].A);
                }//end loop on inner regions
            } //end of loop over regions
	        cout << "total number of hexes counted " << countHex << endl;
            int idx = 0;
            for (auto &h : S[j].Hgrid->hexen) {
                array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionA[idx]);
                iidisp.drawHex(h.position(),(h.d/2.0f),colour);
            }
            cout << "just before redraw display 0" << endl;
            usleep (1000000);
            disp.redrawDisplay();
            usleep (1000000); // one hundred seconds
            disp.saveImage(logpath + "/nnField" + iter + ".png");
            disp.closeDisplay();
        } // end of print on numprint
#endif
    } //end of numsteps loop
//cout << " just after time step i = " << i << endl;
//code run at end of timestepping
//first save the  ofstream outFile;
    morph::HdfData data(fname);
    for (unsigned int j=0;j<numpoints;j++) {
        std::string nstr = "a" + to_string(j);
        char * nst = new char[nstr.length()+1];
        //std::copy(nstr.begin(),nstr.end(),nst);
        std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "b" + to_string(j);
        char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
        std::strcpy(ccst,ccstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        fdata.add_contained_vals(ccst,S[j].B);
        fdata.add_contained_vals(nst,S[j].A);
    }


    cout << " just after writing data "  << endl;
// post run analysis
// look at correlation between adjacent edges

    cout << "before correlate_edges" << endl;
    FLT avAbsCorrelation = M.correlate_edges(0,lZero);
    cout << "after correlate_edges zero morph" << endl;
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
        if (M.innerRegion[j]){
            int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            tempArea = M.regArea(j); //area by hexes
            tempPerimeter = M.regPerimeter(j); //perimeter by hexes
            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].A );
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].A);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            //radial degree
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].A);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " << std::endl;
            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
		radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].A);
		newdegreeRadius = L.find_zeroRadius(radiusVector,3);
		if (newdegreeRadius > degreeRadius) {
                    degreeRadius = newdegreeRadius;
		}
            }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


            regionCount++;
        }//end of if on non-zero regions
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
        if (M.innerRegion[j]){
            countRegions++;
            cout << "just before regAFrac" << endl;
            occupancy += M.regnnfrac(j);
            cout << "just after regAFrac" << endl;
            tempArea = M.regArea(j);
            tempPerimeter = M.regPerimeter(j);
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].A);
            cout << "just ager sectorize Dangle" << endl;
            degreeAngle = L.find_zeroDAngle(angleDVector);
            avDegreeAngle += degreeAngle;
            //radial degree
            degreeRadius = 0;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].A);
            cout << "just ager sectorize Dradius" << endl;
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            avDegreeRadius += degreeRadius;
            cout << "just before writing results first time" << endl;
            degfile1 << degreeAngle/2 << " " << degreeRadius << " " << occupancy << " " << tempArea<< " "<< tempPerimeter<< " xval " << M.centres[j].first << " yval " << M.centres[j].second <<endl<<flush;
            cout << "just after writing results first time" << endl;
        }              //end of loop on inner regions
    } //end of loop on NUMPOINTs
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
    occupancy = occupancy / (1.0 * countRegions);
    jfile <<D_A<<" "<<D_B<<" "<<k1<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation<<endl;
//end of integration on the original tesselation


// section for solving on the circle Tessllation
// if (skipMorph) return 0;
    cout << "just before creating circle Tessellation" <<endl;
    vector<ksSolver> S;
    // do we perturb the centres?
    vector<FLT> rad;
    rad.resize(0);
    vector<FLT> maxRadius;
    vector<FLT> minRadius;
    //set the vector for perturbation of the centres
    maxRadius.resize(0);
    minRadius.resize(0);
    for (int j=0; j<NUMPOINTS; j++) {
        maxRadius.push_back(M.max_radius(j,true));
        minRadius.push_back(M.min_radius(j,true));
    }
    cout << "just before perturbing centres 0 " <<endl;
    if (Lperturb) {
        for (int j=0; j<NUMPOINTS; j++) {
            rad.push_back(minRadius[j]*3.0*ruf.get());
        }
    }
    cout << "just after perturbing centres 0 " <<endl;

    vector<FLT> D_AVal;
    D_AVal.resize(NUMPOINTS,0.0);
    vector<FLT> D_BVal;
    D_BVal.resize(NUMPOINTS,0.0);
    vector<FLT> k1Val;
    k1Val.resize(NUMPOINTS,0.0);
    vector<FLT> k2Val;
    k2Val.resize(NUMPOINTS,0.0);
    vector<FLT> k3Val;
    k3Val.resize(NUMPOINTS,0.0);
    vector<FLT> k4Val;
    k4Val.resize(NUMPOINTS,0.0);
    vector<pair<FLT,FLT>> perturbCentres;
    std::pair<FLT, FLT> centroid;
    if (Lperturb) {
        for (int j = 0;j<NUMPOINTS;j++) {
            FLT angle = 2.0 * PI * ruf.get();
            centroid =  M.baryCentre(j);
            cout << "just before perturbing centres 1" <<endl;
            centroid.first += rad[j]*cos(angle);
            centroid.second += rad[j]*sin(angle);
            cout << "just before perturbing centres 2" <<endl;
            perturbCentres.push_back(centroid);
        }
        cout << "just after perturbing centres 1" <<endl;
    }
    else {
        for (int j=0; j<NUMPOINTS; j++) {
            centroid =  M.baryCentre(j);
            perturbCentres.push_back(centroid);
        }
    }
    cout << "just after perturbing centres 2" <<endl;
    for (int j=0; j<NUMPOINTS; j++) {
        FLT radius = 0;
        if (LDn) {
            radius = maxRadius[j] + 3.0 * minRadius[j];
            FLT area = M.hexArea*M.regArea(j);
            D_AVal[j] = D_A * sqrt((PI * radius * radius / area));
            D_BVal[j] = D_B * sqrt((PI * radius * radius / area));
            k1Val[j] = k1 * sqrt((PI * radius * radius / area));
            k2Val[j] = k2 * sqrt((PI * radius * radius / area));
            k3Val[j] = k3 * sqrt((PI * radius * radius / area));
            k4Val[j] = k4 * sqrt((PI * radius * radius / area));
            cout << "D_AVal region " << j << " = " << D_AVal[j] << " D_A = " << D_A << " PI " << PI << endl;
        }
        else {
            radius = minRadius[j];
            D_AVal[j] = D_A;
            D_BVal[j] = D_B;
            k1Val[j] = k1;
            k2Val[j] = k2;
            k3Val[j] = k3;
            k4Val[j] = k4;
            cout << "D_AVal region " << j << " = " << D_AVal[j] << " D_A = " << D_A << " PI " << PI << endl;
        }

        cout << " radius = " << radius << endl;
        S.push_back(ksSolver(scale, xspan, logpath, radius, perturbCentres[j], lengthScale));
        cout << "in the loop populating the ksVector"<< j << " centroid.x " << perturbCentres[j].first - M.baryCentre(j).first << " centroid.y " << perturbCentres[j].second - M.baryCentre(j).second << endl;

    } //end of loop over NUMPOINTs setting D_A etc..
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

#ifdef COMPILE_PLOTTING
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
#endif
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
            std::string ccstr = "b" + to_string(j);
            cout << " j string " << to_string(j) << " length" << ccstr.length()<< endl;
            char * ccst = new char[ccstr.length()+1];
            std::strcpy(ccst,ccstr.c_str());
            std::string nstr = "a" + to_string(j);
            char * nst = new char[nstr.length()+1];
            std::strcpy(nst,nstr.c_str());
            cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
            ginput.read_contained_vals(ccst,S[j].B);
            ginput.read_contained_vals(nst,S[j].A);
        }
    }
    else
    {
        for (unsigned int j=0;j<NUMPOINTS;j++) {
            morph::RandUniform<FLT> ruf1(seed);
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                FLT choice = ruf.get();
                if (choice > 0.5)   {
                    S[j].A[h.vi] = - ruf1.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf1.get() * aNoiseGain + ccInitialOffset;
                }
		else	{
                    S[j].A[h.vi] =  ruf1.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] =  ruf1.get() * aNoiseGain + ccInitialOffset;
		}
                /*
                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                    S[j].A[h.vi] = (S[j].A[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                    S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
                } //end of if on boundary distance
                */
            }//end of loop over region
            usleep(1000);
        }//end of loop over all regions
    } //end of else on Lcontinue
    // begin time stepping loop cutter grid solved via ksSolver
    // set up vectors for determining convergence

    std::vector<FLT> Adiff;
    std::vector<std::vector<FLT>> Apre;
    std::vector<std::vector<FLT>> Acurr;
    Adiff.resize(NUMPOINTS);
    Apre.resize(NUMPOINTS);
    Acurr.resize(NUMPOINTS);
    FLT AdiffSum = 0.0;

    //initilise all Apre vectors to above possible field
    for (unsigned int j=0; j<NUMPOINTS;j++) {
        Apre[j].resize(S[j].A.size(),1000.0);
    }
    //now set the parameters for each solver
    for (int j=0; j<numpoints;j++) {
        S[j].setParams(D_A, D_B, k1, k2, k3, k4);
    }
    //now set the parameters for each solver
    for (int j=0; j<numpoints;j++) {
        S[j].setParams(D_AVal[j], D_BVal[j], k1Val[j], k2Val[j], k3Val[j], k4Val[j]);
    }
    //start of time-stepping loo
    for (int i=0;i<numsteps;i++) {
        for (int j = 0;j<NUMPOINTS;j++) //loop over regions
        {
            S[j].step(dt);
            if (i%checkEvery == 0) {
                AdiffSum = 0.0;
                Acurr[j] = S[j].A;
                Adiff[j] = L.normedDiff(Apre[j], Acurr[j]);
                cout << "nomm Apre " << L.vectNorm(Apre[j]) << " normCurr " << L.vectNorm(Acurr[j]) << " diff " << Adiff[j] << endl;
                Apre[j] = Acurr[j];
                AdiffSum += Adiff[j];
            }
        }
        if (AdiffSum/(NUMPOINTS*1.0) < diffTol) {
            cout << "morphed converged step " << i << " field diff " << AdiffSum/(1.0*NUMPOINTS) << " diffTol " << diffTol << std::endl;
            break;
        }
    }
    //now populate regions with the final field values
    for (int j=0; j<NUMPOINTS; j++) {
        M.A[j] = S[j].A;
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
        vector<FLT> normalA = M.meanzero_vector(M.sortedBoundaryA[j]);
        //M.printFLTVect(fileString1,normalA);
        cout << " sortedBundary size " << M.sortedBoundary[j].size() << " normalA size " << normalA.size() << endl;
        zeroIndices = L.find_zeroIndices(normalA);
        cout << "zeroIndices size" << zeroIndices.size() << endl;
        int zsize = zeroIndices.size();
        for (int i=0; i<zsize; i++) {
            zeroPhi.push_back(M.sortedBoundaryPhi[j][zeroIndices[i]]);
        }
        M.printIntVect(fileString, zeroIndices);
        M.printFLTVect(fileString1, zeroPhi);
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
#ifdef COMPILE_PLOTTING
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
            for (auto h : S[j].Hgrid->hexen) {
                if (h.testFlags(HEX_IS_REGION_BOUNDARY) == true) {
                    ndisp.drawHex (h.position(), (h.d/2.0f), cl_d);
                    boundaryCount++;
                }
            }
            cout << " region " << j << " has boundary length " << boundaryCount << endl;
            } //end of loop on inner regions
        }//end of loop on NUMPOINTS to plot boundaries
        usleep (1000000);
        cout << "before redrawDisplay 2 " << endl;
        ndisp.redrawDisplay();
        cout << "after redrawDisplay 2" << endl;
        usleep (1000000); // one hundred seconds
        ndisp.saveImage(logpath + "/Tessellation2" + iter + ".png");
        ndisp.closeDisplay();
        //set up display
        morph::Gdisplay mmdisp(window, window, 0, 0, "Laplacian Morphed", rhoInit, 0.0, 0.0);
        mmdisp.resetDisplay (fix, eye, rot);
        for (int j = 0;j<NUMPOINTS;j++) {
            if (M.innerRegion[j]) {
                cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
                cout << " region " << j << " size is " << regsize << endl;
                vector<FLT> tempA;
                vector<FLT> regionA;
        	for (auto h : S[j].Hgrid->hexen) {
                    if (M.hexInRegion(j,h)) {
                        tempA.push_back(S[j].lapA[h.vi]);
                        cout << "tempA " << S[j].lapA[h.vi] << " h.vi " << h.vi  << endl;
                    }
        	}
        //normalise over the region then write normalised values to normalised A over hexGrid
                regionA = L.normalise(tempA);
                int idx = 0;
        	for (auto h : S[j].Hgrid->hexen)
        	{
                    if (M.hexInRegion(j,h)) {
   			array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionA[idx]);
                        mmdisp.drawHex(h.position(),(h.d/2.0f),colour);
                        cout << "just after drawHex 2  value " << regionA[h.vi] << endl;
                        /*
                        std::array<FLT,3> pos2 = {perturbCentres[j].first,perturbCentres[j].second,0};
                        mmdisp.drawHex (pos2,offset,4.0*hexWidth,cl_d);
                        */
                        idx++;
                    }
                }
                cout << "idx in new Hex routine  " << idx << endl;
            } //end of loop on inner regions
        } //end of loop over regions
        cout << "just before redraw display 1" << endl;
        mmdisp.redrawDisplay();
        cout << "just after redraw display 1" << endl;
        cout << "just after to_string"<<endl;
        mmdisp.saveImage(logpath + "/nnLapl2" + iter + ".png");
        usleep (1000000); // one hundred seconds
        cout << "just after saveImage 1" << endl;
        mmdisp.closeDisplay();
    //set up display
        morph::Gdisplay nndisp(window, window, 0, 0, "A field", rhoInit, 0.0, 0.0);
        nndisp.resetDisplay (fix, eye, rot);
        for (int j = 0;j<NUMPOINTS;j++) {
            if (M.innerRegion[j]) {
                cout << "in print routine"<<endl;
                unsigned int regsize = S[j].Hgrid->num();
                cout << " region " << j << " size is " << regsize << endl;
                vector<FLT> tempA;
                vector<FLT> regionA;
        	for (auto h : S[j].Hgrid->hexen) {
                    if (M.hexInRegion(j,h)) {
                        tempA.push_back(S[j].A[h.vi]);
                        cout << "tempA " << S[j].A[h.vi] << " h.vi " << h.vi  << endl;
                    }
        	}
        //normalise over the region then write normalised values to normalised A over hexGrid
                regionA = L.normalise(tempA);
                int idx = 0;
        	for (auto h : S[j].Hgrid->hexen)
        	{
                    if (M.hexInRegion(j,h)) {
   			array<FLT,3> colour = morph::Gdisplay::getJetColorF(regionA[idx]);
                        nndisp.drawHex(h.position(),(h.d/2.0f),colour);
                        cout << "just after drawHex 2  value " << regionA[h.vi] << endl;
                        /*
                        std::array<FLT,3> pos2 = {perturbCentres[j].first,perturbCentres[j].second,0};
                        mmdisp.drawHex (pos2,offset,4.0*hexWidth,cl_d);
                        */
                        idx++;
                    }
                }
                cout << "idx in new Hex routine  " << idx << endl;
            } //end of loop on inner regions
        } //end of loop over regions
        cout << "just before redraw display 1" << endl;
        nndisp.redrawDisplay();
        cout << "just after redraw display 1" << endl;
        cout << "just after to_string"<<endl;
        nndisp.saveImage(logpath + "/nnField2" + iter + ".png");
        usleep (1000000); // one hundred seconds
        cout << "just after saveImage 1" << endl;
        nndisp.closeDisplay();
    } //end of printing of gingerbread cut field
#endif
    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData gdata(gname);
    for (unsigned int j=0;j<NUMPOINTS;j++) {
	std::string nstr = "a" + to_string(j);
	char * nst = new char[nstr.length()+1];
	//std::copy(nstr.begin(),nstr.end(),nst);
	std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "b" + to_string(j);
        char * ccst = new char[ccstr.length()+1];
        std::strcpy(ccst,ccstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        gdata.add_contained_vals(ccst,S[j].B);
        gdata.add_contained_vals(nst,S[j].A);
    }
    cout << " just after writing data "  << endl;

    gfile << endl << "analysis on first morphing iteration " << endl;
    for (int j=0;j<NUMPOINTS;j++) {
        if (M.innerRegion[j]){
            int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);
            // digital version
            angleDVector = M.sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors,S[j].A);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            // analogue version
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors,S[j].A);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M.find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            //radial degree
            degreeRadius = -100;
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].A);
            //gfile <<"after sectorize_reg_radius"<<endl;
            // radiusVector = M.meanzero_vector(radiusVector);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;
            ///radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
                radiusVector = M.sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].A);
                newdegreeRadius = L.find_zeroRadius(radiusVector,3);
                if (newdegreeRadius > degreeRadius)
                    degreeRadius = newdegreeRadius;
                }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
            degfile1 << degreeAngle/2 << " " << degreeRadius << " " << occupancy << " " << tempArea<< " "<< tempPerimeter<< " xval " << M.centres[j].first << " yval " << M.centres[j].second <<endl<<flush;
            cout << "just after writing results first time" << endl;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINT

    if (skipMorph) return 0;
    // clear the edges
    M.edges_clear();
    // redissect the boundaries
    cout << "just before renewDissect morph2 " << endl;
    for (int j=0;j<NUMPOINTS;j++) {
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
//    avAbsCorrelation = 0;
    avAbsCorrelation = M.correlate_edges(1);
    cout << "just after M.correlate_edges, third pass" << endl;
    for (int j=0;j<NUMPOINTS;j++) {
        if (M.innerRegion[j]) {
	    countRegions++;
            //avAbsCorrelation += M.renewcorrelate_edges(j,1,lZero);
	} //end of loop on NUMPOINTs
    }
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
    cout << D_A <<"  "<< D_B <<" "<< k1 <<" " << avAbsCorrelation << " after morphing  2" << endl;
    // write the edge vectors all interpolated to a uniform size
    std::map<int, vector<FLT>>::iterator ptr;
    std::string sideA = logpath + "/edgeA.data";
    int ecount = 0;
    for (ptr = M.edgeA.begin(); ptr != M.edgeA.end(); ptr++) {
        vector<FLT> tempVect;
        tempVect = M.normaliseVect(ptr->second);
        M.printFLTVect(sideA, tempVect);
        cout << "edgeA key " << ptr->first << " element " << ecount << " vector size " << ptr->second.size() << " edges size " << M.edges[ptr->first].size() <<  endl;
        ecount++;
    }
//end of integration after solving on polygonal regions via ksSolver
    return 0;
}
