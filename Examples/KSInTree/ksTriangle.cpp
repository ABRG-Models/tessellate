/*
using morph::HexGrid;
using morph::HdfData;:
using morph::Tools;
*/
#include "ksRegion.h"
#include "analysis.h"
#include <morph/rngd.h>
using namespace std;



#ifdef COMPILE_PLOTTING
/*!
 * If COMPILE_PLOTTING is defined at compile time, then include the display and
 * plotting code. I usually put all the plotting code inside #defines like this so
 * that I can compile a version of the binary without plotting, for parameter searches
 * in which I am only going to be saving out HDF5 data.
 */
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <morph/ColourMap.h>
#include <morph/VisualDataModel.h>
#include <morph/Scale.h>
#include <morph/Vector.h>

//! Helper function to save PNG images with a suitable name
void savePngs (const std::string& logpath, const std::string& name,
               unsigned int frameN, morph::Visual& v)
{
    std::stringstream ff1;
    ff1 << logpath << "/" << name<< "_";
    ff1 << std::setw(5) << std::setfill('0') << frameN;
    ff1 << ".png";
    v.saveImage (ff1.str());
}
#endif

# include <morph/Config.h>
/*
 * using directives just before main()
 */
#ifdef COMPILE_PLOTTING
using morph::Visual;
using morph::HexGridVisual;
using morph::ColourMap;
using morph::ColourMapType;
using morph::VisualDataModel;
using morph::Scale;
using morph::Vector;
#endif
using morph::Config;
using morph::Tools;
using std::string;
using std::stringstream;
using std::cerr;
using std::cout;
using std::endl;
using std::vector;
using std::chrono::milliseconds;
using std::chrono::duration_cast;
using std::chrono::steady_clock;

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    string jsonfile = argv[1];
    vector<morph::BezCurvePath<float>> triangleBound;
    triangleBound.resize(0);
    hexGeometry* hGeom;
    double dTress = 0.5;
    // open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
    cout << "just before reading parameters" << endl;
    double dt = conf.getDouble("dt",0.0002);
    cout << "dt " << dt << endl;
    double Dn = conf.getDouble("Dn",36.0);
    cout << "Dn " << Dn << endl;
    double Dchi = conf.getDouble("DChi",0.0);
    cout << "Dchi " << Dchi << endl;
    double Dc = conf.getDouble("Dc",12.0);
    cout << "Dc  " << Dc << endl;
//    double ds = conf.getInt("ds",0.01);
//    cout << "ds " << ds << endl;
    int scale = conf.getInt("scale", 8);
    cout << "scale " << scale << endl;
    double xspan = conf.getDouble("xspan",4.0);
    cout << "xspan " << xspan << endl;
    int numsteps = conf.getInt("numsteps",200);
    cout << "numsteps " << numsteps << endl;
    int numAdjust = conf.getInt("numAdjust",101);
    cout << "numAdjust " << numAdjust << endl;
    int numprint = conf.getInt("numprint",85);
    cout << "numprint " << numprint << endl;
    string logpath = conf.getString("logpath","../logs");
    cout <<  logpath << endl;
    bool Lgraphics = conf.getBool("Lgraphics", true);
    cout << "Lgraphics " << Lgraphics << endl;
    bool LDn = conf.getBool("LDn", true);
    cout << "LDn " << LDn << endl;
    double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
    cout << "boundarFalloffDist " << boundaryFalloffDist << endl;
    double aNoiseGain = conf.getDouble("aNoiseGain",0.2);
    cout << "aNoiseGain " << aNoiseGain << endl;
    int numSectors = conf.getInt("numsectors",14);
    cout << "numSectors " << numSectors << endl;
    bool Lcontinue = conf.getBool("Lcontinue",false);
    cout << "Lcontinue " << Lcontinue << endl;
    double nnInitialOffset = conf.getDouble("nnInitialOffset", 1.0);
    cout << "nnInitialOffset " << nnInitialOffset << endl;
    double ccInitialOffset = conf.getDouble("ccInitialOffset",2.5);
    cout << "ccInitialOffset " << ccInitialOffset << endl;
    bool overwrite_logs = conf.getBool("overwrite_logs",true);
    cout << "overwrite_logs " << overwrite_logs << endl;
    bool skipMorph  = conf.getBool("skipMorph",true);
    cout << "skipMorph " << skipMorph << endl;
    bool lPerturb = conf.getBool("lPerturb",true);
    cout << "lPerturb " << lPerturb << endl;
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    cout << "LfixedSeed " << LfixedSeed << endl;
    double ratio = conf.getDouble("ratio", 1.0);
    cout << "ratio  " << ratio  << endl;
    double pRatio = conf.getDouble("pRatio",0.0);
    cout << "pRatio " << pRatio << endl;
    bool Lrows = conf.getBool("Lrows", true);
    cout << "Lrows " << Lrows << endl;
    bool Lzero = conf.getBool("Lzero",true);
    cout << "Lzero " << Lzero << endl;
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << " pRatio " << pRatio << " ratio " << ratio << "  " << logpath <<  endl;
    ofstream afile (logpath + "/centroids.out",ios::app);
    int numpoints;

#ifdef COMPILE_PLOTTING
    // Parameters from the config that apply only to plotting:
  //  const unsigned int plotevery = conf.getUInt ("plotevery", 10);
    // Should the plots be saved as png images?
    const bool saveplots = conf.getBool ("saveplots", false);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    const bool vidframes = conf.getBool ("vidframes", false);
    unsigned int framecount = 0;

    // Window width and height
    const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
    unsigned int win_height_default = static_cast<unsigned int>(0.8824f * (float)win_width);
    const unsigned int win_height = conf.getUInt ("win_height", win_height_default);

    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    Visual v1 (win_width, win_height, "ksTriangle");
    // Set a dark blue background (black is the default). This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v1.bgcolour = {1.0f, 1.0f, 1.0f, 1.0f};

    //set orthographic projection
    v1.ptype = morph::perspective_type::orthographic;
    v1.ortho_bl = {-0.9f, -0.9f};
    v1.ortho_tr = {0.9f, 0.9f};
    // You can tweak the near and far clipping planes
    v1.zNear = 0.001;
    v1.zFar = 20;
    // And the field of view of the visual scene.
    v1.fov = 45;
    // You can lock movement of the scene
    v1.sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v1.setZDefault (conf.getFloat ("z_default", -5.0f));
    v1.setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v1.scenetrans_stepsize = 0.5;
    steady_clock::time_point lastrender = steady_clock::now();
#endif

    /*
     * decide if we are using a fixed or
     * run-time dependent seed. Use milliseconds to
     * allow parallel launching
     */

    unsigned int seed;
    chrono::milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    if (LfixedSeed)
        seed = 1;
    else
        seed = static_cast<unsigned int> (ms1.count());

    morph::RandUniform<double> ruf(seed);

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
    ifstream vfile ("./triangle.inp");
    /*
     * analysis tools
     */
    Analysis L;

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate ksSolver
     */
    vector<ksSolver> S;
    vector<pair<double,double>> centres;
    morph::BezCurvePath<float> outer;
    const int rowX = 2; const int rowY = 2;
    //get numpoints from the above row sizes
    numpoints = 2 * (2*rowX + 1) * (2*rowY + 1);
    //next lines for creating a tessellation of generic triangles.
    //if ratio is set to 1.0 we get equilateral triangles
    //if pRatio is > 0.0 and lPerturb == true then we get a tessellation of
    //perturbed triangles.
    vector<vector<hexGeometry::point>> vtxs;
    vtxs.resize(numpoints);
    vector<hexGeometry::point> vertices;
    Tessellation* M;
    centres.resize(0);
    vector<vector<int>> nbrList;
    if (Lrows) {
        vertices = hGeom->isosVertices(ratio, rowX, rowY, pRatio, lPerturb);
        vector<vector<int>> vIndices;
        nbrList =hGeom->triangleNeighbors(rowX, rowY, vIndices);
        triangleBound = hGeom->genTriangleTess(rowX, rowY, vertices, vIndices, outer);
        for (int j=0; j<numpoints; j++) {
            vtxs[j].push_back(vertices[vIndices[j][0]]);
            vtxs[j].push_back(vertices[vIndices[j][1]]);
            vtxs[j].push_back(vertices[vIndices[j][2]]);
        }
        vector<vector<hexGeometry::point>>::iterator itr;
        cout<< "before fill centres size of vtxs " << vtxs.size() << endl;
        //create the centres of the triangles
        for (itr = vtxs.begin(); itr != vtxs.end(); itr++) {
            cout << " at head of vtxs loop " << endl;
            cout << "size of vtxs first vector " << (*itr).size() << endl;
            cout << (*itr)[0].first << " , " << (*itr)[1].first << " , " << (*itr)[2].first << endl;
            double px = ((*itr)[0].first + (*itr)[1].first + (*itr)[2].first) / 3.0;
            double py = ((*itr)[0].second + (*itr)[1].second + (*itr)[2].second) / 3.0;
            cout << "before centres.push_back()" << endl;
            centres.push_back(std::make_pair(px, py));
        }
        M = new Tessellation(scale, xspan, logpath, outer, vtxs, centres, numpoints, Lrows);
    }
    else {
        //next line for creating a tessellation of equilateral triangles
        //I am being lazy heres and hard coding the number of triangles
        //but is OK because I am creating an ensemble of tessellations
        numpoints = 42; //apologies to Douglas Adams
        triangleBound = hGeom->eqTriangleTess(dTress, centres, outer, vtxs, nbrList);
        M = new Tessellation(scale, xspan, logpath, outer, vtxs, centres, numpoints, Lrows);
        //M= new Tessellation (scale, xspan, logpath, outer, centres, numpoints, Lrows);
        //M->vCoords = vtxs;
        cout << "size of triangleBound " << triangleBound.size() << endl;
    }
    cout << "after Tessellation" << endl;

    M->setRegionList(nbrList);
    //M->dissectBoundary();
    cout << "just before setRadialSegments" << endl;
    M->setRadialSegments();
    cout << "after Tessellation size of M->vCoords" << M->vCoords.size() << endl;
    /*
     * Solve the KS equations over the regions
     */
    for (int j = 0;j<numpoints;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, triangleBound[j], centres[j]));
        //S.push_back(ksSolver(scale, xspan, logpath, M->curvedBoundary[j], centroid));
        cout << "in the loop populating the ksVector"<< j <<endl;
    }

// initialise the fields
    string fname = logpath + "/first.h5";
    cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++)
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
        for (unsigned int j=0;j<numpoints;j++)
		{
	        for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		        double choice = ruf.get();
		        if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
				}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

#ifdef COMPILE_PLOTTING

    // Spatial offset, for positioning of HexGridVisuals
    Vector<float> spatOff;
    float xzero = 0.0f;

    // A. Offset in x direction to the left.
    xzero = -0.11;
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    Scale<FLT,float> zscale; zscale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    Scale<FLT,float> cscale; cscale.do_autoscale = true;
    unsigned int Agrid[numpoints];
    for (int j=0; j<numpoints; j++) {
        vector<double> regionNN;
//normalise over the region then write normalised values to normalised NN over hexGrid
        regionNN = L.normalise(S[j].NN);
        Agrid[j] = v1.addVisualModel (new HexGridVisual<FLT> (v1.shaderprog,
                                                                    v1.tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::Jet));
      }

#endif

    int boundaryCount = 0;
    // begin first time stepping loop
    cout << "before time steps dt " << dt << " Dn" << Dn << " Dchi " << Dchi << " Dc " << Dc << endl;
    for (int i=0;i<numsteps;i++)
    {
   	    for (int j = 0;j<numpoints;j++) //loop over regions
        {
     	    S[j].step(dt, Dn, Dchi, Dc);
	    }
#ifdef COMPILE_PLOTTING
	    if (i%numprint == 0)
	    {
          //set up display
   	        for (int j = 0;j<numpoints;j++) //loop over regions
	        {
	 	        cout << "in print routine i "<< i  << " j " << j << endl;
		        vector<double> regionNN;
//normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(S[j].NN);
                VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v1.getVisualModel (Agrid[j]);
                avm->updateData (&regionNN);
                avm->clearAutoscaleColour();
            }
            // rendering the graphics. After each simulation step, check if enough time
            // has elapsed for it to be necessary to call v1.render().
            cout << " before rendering the graphics " << endl;
            steady_clock::duration sincerender = steady_clock::now() - lastrender;
            if (duration_cast<milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
                glfwPollEvents();
                cout << " before rendering the graphics " << endl;
                v1.render();
                cout << " after rendering the graphics " << endl;
                lastrender = steady_clock::now();
            }
            cout << " before vidframes " << endl;
            if (vidframes) {
                savePngs (logpath, "nnField", framecount, v1);
                ++framecount;
            } else {
                  savePngs (logpath, "nnField", 0, v1);
            }
        }//end of loop on numprint drawing fields
#endif
    } // end of first time stepping
    //v1.keepOpen();

//code run at end of timestepping

    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData fdata(fname);
	for (unsigned int j=0;j<numpoints;j++)
	{
		std::string nstr = "n" + to_string(j);
	    char * nst = new char[nstr.length()+1];
		//std::copy(nstr.begin(),nstr.end(),nst);
		std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
		std::strcpy(ccst,ccstr.c_str());
		cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        fdata.add_contained_vals(ccst,S[j].CC);
        fdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M->X[0]);
        //data.add_contained_vals("Y",M->X[1]);
    }
    fdata.add_val ("/Dchi", Dchi);
    fdata.add_val ("/Dn", Dn);
    fdata.add_val ("/Dc",Dc);

    /*
     * Now populate the Tessellation regions with NN
     */
    cout << "first setting of centroids" << endl;
    // repopulate the regions
    for (int j=0;j<numpoints;j++)
    {
        double NNNorm;
        cout << "jsut before renewRegion J " << j << endl;
        NNNorm = M->renewRegion(j,S[j].Hgrid->hexen, S[j].NN);
        cout << "nnNorm after M->renewRegion " << NNNorm << endl;
    }
    for (int j=0;j<numpoints;j++)
    {
        M->renewBoundary(j,S[j].Hgrid->hexen);
    }
    cout << "just before edges_clear() first time" << endl;
    M->edges_clear();
    cout << "just before setting regionList 1" << endl;
    cout << "just before setting regionList 2" << endl;
    //M->setregionList(nbrList);
    cout << "just before renewDissect " << endl;
    for (int j=0;j<numpoints;j++)
    {
        M->renewDissect(j,0);
    }
    cout << "Edges size " << M->edges.size() << endl;
    int xcount=0;
    for (auto itr = M->edges.begin(); itr != M->edges.end(); itr++) {
         xcount++;
         cout << "edge " << xcount << " key " << itr->first << " list " << endl;
         M->printIntVect(logpath + "/renewDissect.out", itr->second);
    }
    cout << " renewDissect size of edges " << M->edges.size() << " size of edgeIndex " << M->edgeIndex.size() << endl;


    afile << "centroids for zero morphin " << endl;
    for (int j=0; j<numpoints;j++){
        afile << " Region size " << M->regionHex[j].size() << endl;
    }
    vector <int> radiusDVector;
    vector <int> angleDVector;
	vector <double> angleVector;
    vector <double> radiusVector;
    int degreeRadius;
    int degreeAngle;
	double tempArea = 0;
	double tempPerimeter = 0;
    int angleOffset = 0;
    int radiusOffset = 0;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    degreeRadius = 0;
    degreeAngle = 0;
	double avDegreeAngle = 0;
	double avDegreeRadius = 0;
    double avAbsCorrelation = 0;
    double occupancy = 0;
    int countRegions = 0;
    int max_comp = numpoints * 3;

    for (int j=0;j<numpoints;j++) {
        if (M->regArea(j) != 0){
			int regionCount = 0;
            gfile<<"in the degree loop" << endl;
            //angle degree
            tempArea = M->regArea(j);
            tempPerimeter = M->renewRegPerimeter(j);
            gfile << "before sectorize Dangle" << endl;
            // digital version
            angleDVector = M->sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            // analogue version
            angleVector = M->sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M->find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

            //radial degree
            degreeRadius = -100;
            radiusDVector = M->sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            // radiusVector = M->meanzero_vector(radiusVector);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
            for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
			    radiusVector = M->sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			    newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			    if (newdegreeRadius > degreeRadius) {
                    degreeRadius = newdegreeRadius;
                }
		    }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
        } //end of if on non-zero regions
    } //end of loop on NUMPOINtS


    avDegreeAngle = 0;
    avDegreeRadius = 0;
    occupancy = 0;
    countRegions = 0;
    tempArea = 0;
    tempPerimeter = 0;
    avAbsCorrelation = 0;
    cout << "just before correlate_edges morph1 " << endl;
    M->correlate_edges();
    cout << "just after  correlate_edges morph1 " << endl;
    cout << "just before randomcorrelate_edges morph1 " << endl;
    M->random_correlate(max_comp,0,Lzero);
    for (int j=0;j<numpoints;j++) {
        countRegions++;
        occupancy += M->regNNfrac(j);
        cout << "just after regNNfrac " << endl;
        //tempArea = M->regArea(j);
        //tempPerimeter = M->regPerimeter(j);
        cout << "just before renewcorrelate_edges morph1 " << endl;
        avAbsCorrelation += M->renewcorrelate_edges(j,0,Lzero);
        cout << "just after renewcorrelate_edges morph1 " << endl;
        angleDVector = M->sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
        degreeAngle = L.find_zeroDAngle(angleDVector);
        avDegreeAngle += degreeAngle;
        //radial degree
        degreeRadius = 0;
        radiusDVector = M->sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
        degreeRadius = L.find_zeroDRadius(radiusDVector);
        avDegreeRadius += degreeRadius;

        degfile1 << degreeAngle/2 << " " << degreeRadius << " " << M->regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
    }//end of loop on numpoints
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
    occupancy = occupancy / (1.0 * countRegions);
    cout << "Regions counted in first analysis loop " << countRegions << endl;
    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
    //end of integration after solving on polygonal regions via ksSolver
    //
    //
    if(skipMorph) return 0;
    //Now solve on an array of circles that are circumscribed on the triangles
    S.clear();
#ifdef COMPILE_PLOTTING
   // Before starting the simulation, create the HexGridVisuals.
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    Visual v2 (win_width, win_height, "ksTriangle");
    // Set a dark blue background (black is the default). This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v2.bgcolour = {1.0f, 1.0f, 1.0f, 1.0f};
    // You can tweak the near and far clipping planes
    v2.zNear = 0.001;
    v2.zFar = 20;
    // And the field of view of the visual scene.
    v2.fov = 45;
    // You can lock movement of the scene
    v2.sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v2.setZDefault (conf.getFloat ("z_default", -5.0f));
    v2.setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v2.scenetrans_stepsize = 0.5;
    // if using plotting, then set up the render clock
    lastrender = steady_clock::now();
#endif

    double radius = 1.0 / (1.73 * (2 * rowX + 1));
    vector<double> rad;
    rad.resize(0);
    for (int j=0; j<numpoints; j++) {
        rad.push_back(radius*ruf.get()*0.5);
    }
    cout << " radius = " << radius << endl;
    radius = 2.0*radius;
    for (int j=0;j<numpoints;j++) {
        std::pair<double,double> regionCentre;
        double temp = ruf.get();
        if ((temp >= 0.0) && (temp < 0.25)) {
            regionCentre.first = M->centres[j].first + rad[j];
            regionCentre.second = M->centres[j].second + rad[j];
        }
        if ((temp >= 0.25) && (temp < 0.5)) {
            regionCentre.first = M->centres[j].first + rad[j];
            regionCentre.second = M->centres[j].second - rad[j];
        }
        if ((temp >= 0.5) && (temp < 0.75)) {
            regionCentre.first = M->centres[j].first - rad[j];
            regionCentre.second = M->centres[j].second + rad[j];
        }
        if ((temp >= 0.75) && (temp <= 0.25)) {
            regionCentre.first = M->centres[j].first - rad[j];
            regionCentre.second = M->centres[j].second - rad[j];
        }
        S.push_back(ksSolver(scale, xspan, logpath, radius, regionCentre));
        cout << "in the loop populating the ksVector"<< j << " centroid.x " <<  M->centres[j].first << " centroid.y " <<  M->centres[j].second << endl;

    }
	cout << "just after populating the ksVector"<<endl;
    // repopulate the regions and their boundaries with hex indices
    for (int j=0;j<numpoints;j++)
    {
        M->renewRegion(j,S[j].Hgrid->hexen);
        int rsize = M->regionHex[j].size();
        int ssize = M->regionBound[j].size();
        std::pair bcent = M->baryCentre(j);
        cout << "barycentre of " << j << " x " << bcent.first << " y " << bcent.second << " region size " << rsize << " bound size " << ssize << endl;
    }
    for (int j=0;j<numpoints;j++) {
        M->renewBoundary(j,S[j].Hgrid->hexen);
        int rsize = M->regionHex[j].size();
        int ssize = M->regionBound[j].size();
        std::pair bcent = M->baryCentre(j);
        cout << "barycentre of " << j << " x " << bcent.first << " y " << bcent.second << " region size " << rsize << " bound size " << ssize << endl;
     }
// initialise the fields
    string ggname = logpath + "/second.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(ggname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++)
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
        for (unsigned int j=0;j<numpoints;j++)
		{
            unsigned int seed;
            milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
            seed = static_cast<unsigned int> (ms1.count());
            //seed = j + off;
            morph::RandUniform<double> ruf1(seed);
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                double choice = morph::randDouble();
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
            }//end of loop over region
            usleep(1000);
        }//end of loop over all regions
    } //end of else on Lcontinue

// begin second time stepping loop
    cout << "just beore second time steping loop " << endl;
#ifdef COMPILE_PLOTTING
    unsigned int Bgrid[numpoints];
    for (int j=0; j<numpoints; j++) {
        vector<double> regionNN;
//normalise over the region then write normalised values to normalised NN over hexGrid
        regionNN = L.normalise(S[j].NN);
        Bgrid[j] = v2.addVisualModel (new HexGridVisual<FLT> (v2.shaderprog,
                                                                    v2.tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::Jet));
      }

#endif

    // begin first time stepping loop
    cout << "before time steps dt " << dt << " Dn" << Dn << " Dchi " << Dchi << " Dc " << Dc << endl;
    //adjust Dn because we have doubled the radius
    Dn = 2.0 * Dn;
    Dc = 0.3*Dn;
    for (int i=0;i<2*numsteps;i++)
    {
   	    for (int j = 0;j<numpoints;j++) //loop over regions
        {
     	    S[j].step(dt, Dn, Dchi, Dc);
	    }
#ifdef COMPILE_PLOTTING
	    if (i%numprint == 0)
	    {
          //set up display
   	        for (int j = 0;j<numpoints;j++) //loop over regions
	        {
	 	        cout << "in print routine"<<endl;
		        vector<double> regionNN;
//normalise over the region then write normalised values to normalised NN over hexGrid
                regionNN = L.normalise(S[j].NN);
                VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v2.getVisualModel (Bgrid[j]);
                avm->updateData (&regionNN);
                avm->clearAutoscaleColour();
            }
            // rendering the graphics. After each simulation step, check if enough time
            // has elapsed for it to be necessary to call v2.render().
            steady_clock::duration sincerender = steady_clock::now() - lastrender;
            if (duration_cast<milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
                glfwPollEvents();
                v2.render();
                lastrender = steady_clock::now();
            }

            if (vidframes) {
                savePngs (logpath, "nnFieldC", framecount, v2);
                ++framecount;
            } else {
                  savePngs (logpath, "nnFieldC", 0, v2);
            }
        }//end of loop on numprint drawing fields
#endif
    } // end of second time stepping

    //now populate regions with the final field values
    cout << " just before populating regions " << endl;
    for (int j=0; j<numpoints; j++) {
        M->NN[j] = S[j].NN;
    }
    //not the circle
    cout << "just before the cookie cutter routine" << endl;
    vector<vector<hexGeometry::lineSegment>> regionLineSegments;
    regionLineSegments.resize(numpoints);
    for (int j=0;j<numpoints;j++)
    {
        cout << " just before regionLineSegments " << j << " regionSize " <<  regionLineSegments.size()  << endl;
        regionLineSegments[j] = M->polygonSides(j, true);
        unsigned int sideSize = regionLineSegments[j].size();
        cout << "sideSize " << sideSize << " region " << j << endl;
        int count = 0;
        for (auto& h : S[j].Hgrid->hexen) {
            for (unsigned int i=0; i<sideSize; i++) {
                if (M->hGeo->hexIntersectLineSegment(regionLineSegments[j][i],h)){
                    h.setFlags(HEX_IS_REGION_BOUNDARY);
                    M->regionBound[j].push_back(h);
                    count++;
                }
            }
        }
        cout << "in region " << j << " size of region boundary " << count << endl;
    }
    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData ggdata(ggname);
    for (unsigned int j=0;j<numpoints;j++) {
	std::string nstr = "n" + to_string(j);
	char * nst = new char[nstr.length()+1];
	//std::copy(nstr.begin(),nstr.end(),nst);
	std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
        char * ccst = new char[ccstr.length()+1];
        std::strcpy(ccst,ccstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        ggdata.add_contained_vals(ccst,S[j].CC);
        ggdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M->X[0]);
        //data.add_contained_vals("Y",M->X[1]);
     }
     ggdata.add_val ("/Dchi", Dchi);
     ggdata.add_val ("/Dn", Dn);
     ggdata.add_val ("/Dc",Dc);

     cout << " just after writing data "  << endl;

      gfile << endl << "analysis on first morphing iteration " << endl;

    // clear the edges
    M->edges_clear();
    // redissect the boundaries
    cout << "just before renewDissect morph2 " << endl;
    for (int j=0;j<numpoints;j++)
	{
	    M->renewDissect(j, 1);
    }

    cout << "just after renewDissect morph2 " << endl;
//computing average values over all regions
    avAbsCorrelation = 0;
    tempPerimeter = 0;
    tempArea = 0;
    countRegions = 0;
    M->random_correlate(max_comp, 1, Lzero);
    cout << "just after random correlate_edges morph2 " << endl;
    avAbsCorrelation = 0;
    cout << "just after M->correlate_edges, third pass" << endl;
    for (int j=0;j<numpoints;j++)
	{
	    countRegions++;
        avAbsCorrelation += M->renewcorrelate_edges(j,1,Lzero);
	} //end of loop on NUMPOINTs
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	// jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	cout << Dn <<"  "<< Dchi <<" "<< Dc <<" " << avAbsCorrelation << " after CookieCutter correlation " << endl;



    /*
     * now create an array of morphed regions, first morph
     */

    if (skipMorph) return 0;
    cout << "just before setting curved boundaries first morph" <<endl;
// section for solving on the curved boundaries
// now set the boundaries for the regions, stored in curvedBoundary
    M->populateBoundCurve(1);
    cout << "just after setting curved boundaries " << M->curvedBoundary.size()<<endl;
    for (int j = 0;j<numpoints;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M->curvedBoundary[j], M->centroids[j]));
        cout << "in the loop populating the ksVector"<< j <<endl;
    }
    afile << "first morph  setting of centroids" << endl;
    for (int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M->centroids[j].first << " , " << M->centroids[j].second << " )" << endl;
    }

// initialise the fields
    string gname = logpath + "/second.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    if (Lcontinue)
	{
        morph::HdfData ginput(gname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++)
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
        for (unsigned int j=0;j<numpoints;j++)
		{
	        for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
		        double choice = ruf.get();
		        if (choice > 0.5)
				{
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
				}
		        else
				{
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
				}
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
		 } //end of if on boundary distance
	    }//end of loop over region
	   }//end of loop over all regions
      } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

    // begin second time stepping loop
    for (int i=0;i<numsteps;i++)
    {
	    int countHex = 0;
   	    for (int j = 0;j<numpoints;j++) //loop over regions
        {
     	    S[j].step(dt, Dn, Dchi, Dc);
		    if (i%100 == 0)
		    cout << "just after morph 1 octopus i "<< j << " NN[5] " << S[j].NN[5] << endl;
	    }
		 //cout << "iteration " << i << " of morph 1 time step" << endl;
    } // end of second time stepping

    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData gdata(gname);
	for (unsigned int j=0;j<numpoints;j++)
	{
		std::string nstr = "n" + to_string(j);
	    char * nst = new char[nstr.length()+1];
		//std::copy(nstr.begin(),nstr.end(),nst);
		std::strcpy(nst,nstr.c_str());
    	std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
	//	std::copy(ccstr.begin(),ccstr.end(),cst);
		std::strcpy(ccst,ccstr.c_str());
		cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        gdata.add_contained_vals(ccst,S[j].CC);
        gdata.add_contained_vals(nst,S[j].NN);
        //data.add_contained_vals("X",M->X[0]);
        //data.add_contained_vals("Y",M->X[1]);
     }
     gdata.add_val ("/Dchi", Dchi);
     gdata.add_val ("/Dn", Dn);
     gdata.add_val ("/Dc",Dc);
     cout << " just after writing data "  << endl;
     // write the NN and CC vals for each region
      gfile << endl << "analysis on first morphing iteration " << endl;   // repopulate the regions
    for (int j=0;j<numpoints;j++)
    {
        M->renewRegion(j,S[j].Hgrid->hexen,S[j].NN);
    }
    for (int j=0;j<numpoints;j++)
    {
        M->renewBoundary(j,S[j].Hgrid->hexen);
    }
    afile << "second setting of centroids" << endl;
    for (int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M->centroids[j].first << " , " << M->centroids[j].second << " )" << endl;
    }
	cout << "just after populating the ksVector"<<endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
      for (int j=0;j<numpoints;j++) {
              if (M->regArea(j) != 0){
			      int regionCount = 0;
                  gfile<<"in the degree loop" << endl;
                  //angle degree
                  tempArea = M->regArea(j);
                  tempPerimeter = M->renewRegPerimeter(j);

                  // digital version
                  angleDVector = M->sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  degreeAngle = L.find_zeroDAngle(angleDVector);
                  gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

                  // analogue version
                  angleVector = M->sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
                  angleVector = L.meanzero_vector(angleVector);
                  //degreeAngle = M->find_max(angleVector,3);
                  degreeAngle = L.find_zeroAngle(angleVector,3);
                  gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
                  //radial degree
                  degreeRadius = -100;
                  radiusDVector = M->sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                  //gfile <<"after sectorize_reg_radius"<<endl;
                  // radiusVector = M->meanzero_vector(radiusVector);

                  degreeRadius = L.find_zeroDRadius(radiusDVector);
                  gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl ;

                  ///radial degree
                  degreeRadius = -100;
                  int newdegreeRadius = 0;
                  for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3)
				  {
			          radiusVector = M->sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
			          newdegreeRadius = L.find_zeroRadius(radiusVector,3);
			          if (newdegreeRadius > degreeRadius)
				      degreeRadius = newdegreeRadius;
		          }


                  gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;


                  regionCount++;
        } //end of if on non-zero regions
    } //end of loop on numpoints


    //why does this have no effect?
    // redissect the boundaries
    cout << "just before renewDissect second time" << endl;
    for (int j=0;j<numpoints;j++)
    {
        M->renewDissect(j,1);
    }
    cout << "Edges size " << M->edges.size() << endl;


	      avDegreeAngle = 0;
	      avDegreeRadius = 0;
	      occupancy = 0;
	      countRegions = 0;
		  tempArea = 0;
		  tempPerimeter = 0;
		  avAbsCorrelation = 0;
		  cout << "just after renewcorrelate_edges morph1 " << endl;
		  M->random_correlate(max_comp, 1,Lzero);
		  cout << "just after randomcorrelate_edges morph1 " << endl;
          for (int j=0;j<numpoints;j++) {
	        if (M->regArea(j) != 0)
			{
	          countRegions++;
		      occupancy += M->regNNfrac(j);
              tempArea = M->regArea(j);
              tempPerimeter = M->regPerimeter(j);
              avAbsCorrelation += M->renewcorrelate_edges(j,1,Lzero);
              angleDVector = M->sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
              degreeAngle = L.find_zeroDAngle(angleDVector);
		      avDegreeAngle += degreeAngle;
              //radial degree
		      degreeRadius = 0;
              radiusDVector = M->sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
              degreeRadius = L.find_zeroDRadius(radiusDVector);
              avDegreeRadius += degreeRadius;

              degfile2 << degreeAngle/2 << " " << degreeRadius << " " << M->regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
	       } //end of if on non-zero regions
	    } //end of loop on numpoints
        avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
        avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
        avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
	    occupancy = occupancy / (1.0 * countRegions);
	    // jfile << avDegreeAngle <<" "<<avDegreeRadius<<endl;
	    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
//end of integration after first morphing
//begin second morphing
    if (skipMorph) return 0;
// now set the boundaries for the regions, stored in curvedBoundary
    M->populateBoundCurve(0);
    S.resize(0);
    cout << "just after setting curved boundaries second iteration" << M->curvedBoundary.size()<<endl;
    for (int j = 0;j<numpoints;j++)
    {
        std::pair<float, float> centroid =  M->baryCentre(j);
        S.push_back(ksSolver(scale, xspan, logpath, M->curvedBoundary[j], M->centroids[j]));
        cout << "in the loop populating the ksVector"<< j <<endl;
    }
    cout << "just after populating the ksVector"<<endl;
// initialise the fields
        cout<< "just before first data read"<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
        string hname = logpath + "/third.h5";
        if (Lcontinue) {
             morph::HdfData hinput (hname,1);
             for (int j=0;j<numpoints;j++){
                std::string nstr = "n" + to_string(j);
                char * nst = new char[nstr.length()+1];
                std::strcpy(nst,nstr.c_str());
                std::string ccstr = "c" + to_string(j);
                char * ccst = new char[ccstr.length()+1];
                cout << "labels "<< nstr <<" , " << ccstr<<endl;
                std::strcpy(ccst,ccstr.c_str());
                hinput.read_contained_vals(nst,S[j].NN);
                hinput.read_contained_vals(ccst,S[j].CC);
                cout<< "just after input of NN and CC1"<< endl;
       //   input.close();
            }
	 }
     else {
        for (int j=0;j<numpoints;j++) {
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                double choice = ruf.get();
                if (choice > 0.5)
                {
                    S[j].NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
                }
                else {
                    S[j].NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                    S[j].CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
                }

                if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                    double bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                    S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                    S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
                } //end of if on boundary distance
            }//end of loop over region
        }//end of loop over all regions
     } //end of else on Lcontinue
     cout <<  "just after field creation" << endl;


     //begin second time stepping loop
    int loopsteps = 0;
    for (int i=0;i<numsteps;i++)
	{
	  //cout << " head of second time stepping loop i " << i << endl;
        for (int j = 0;j<numpoints;j++) {
     	    S[j].step(dt, Dn, Dchi, Dc);
     	 //  S[j].step(dt, Dn, Dchi, Dc);
	    }
    }
    //code run at end of timestepping
    //first save the  ofstream outFile;
    morph::HdfData hdata(hname);
	cout <<"after creating hdata "<< hname << endl;
	for (int j=0;j<numpoints;j++) {
		std::string nstr = "n" + to_string(j);
	   	char * nst = new char[nstr.length()+1];
		std::strcpy(nst,nstr.c_str());
		std::string ccstr = "c" + to_string(j);
	    char * ccst = new char[ccstr.length()+1];
		std::strcpy(ccst,ccstr.c_str());
		cout <<"in third hdf5 write "<<nstr <<" , "<<ccstr<<endl;
        hdata.add_contained_vals(ccst,S[j].CC);
        hdata.add_contained_vals(nst,S[j].NN);
	}
    //data.add_contained_vals("X",M->X[0]);
    //data.add_contained_vals("Y",M->X[1]);
    // hdata.add_val ("/Dchi", Dchi);
    // hdata.add_val ("/Dn", Dn);
    // hdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;

    gfile << endl << "analysis on second morphing iteration " << endl;
    // repopulate the regions
    for (int j=0;j<numpoints;j++)
    {
        double result = M->renewRegion(j,S[j].Hgrid->hexen, S[j].NN);
        cout << "after renewRegion total NN" << result;
    }
    for (int j=0;j<numpoints;j++)
    {
        M->renewBoundary(j,S[j].Hgrid->hexen);
    }
    cout << "after renewBoundary " << endl;
//    swap the radialAngles to the vCoords
    M->swapRadialSegments(true);
    cout << "after swapRadialSegments " << endl;
    // redissect the boundaries
    for (int j=0;j<numpoints;j++)
    {
        M->renewDissect(j,2);
    }
    cout << "after renewDissect " << endl;
    cout << "Edges size " << M->edges.size() << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (int j=0;j<numpoints;j++) {
	    if (M->regArea(j) != 0) {
	        int regionCount = 0;
            gfile<<"in the degree loop" << endl;
	        //angle degree
            tempArea = M->regArea(j);
            tempPerimeter = M->renewRegPerimeter(j);
            // digital version
            angleDVector = M->sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            degreeAngle = L.find_zeroDAngle(angleDVector);
            gfile << "region "<< j << " degreeDAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;

           // analogue version
            angleVector = M->sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            angleVector = L.meanzero_vector(angleVector);
            //degreeAngle = M->find_max(angleVector,3);
            degreeAngle = L.find_zeroAngle(angleVector,3);
            gfile << "region "<< j << " degreeAngle "<< degreeAngle << "  " << tempArea<< "  "<< tempPerimeter<<endl<<flush;
            //radial degree
            degreeRadius = -100;
            radiusDVector = M->sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            //gfile <<"after sectorize_reg_radius"<<endl;
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            gfile  << "region "<< j << " degreeDRadius "<< degreeRadius << "  " <<endl;
            //radial degree
            degreeRadius = -100;
            int newdegreeRadius = 0;
			for (int angleOffset=0; angleOffset<numSectors -1; angleOffset += 3) {
                radiusVector = M->sectorize_reg_radius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
                newdegreeRadius = L.find_zeroRadius(radiusVector,3);
                if (newdegreeRadius > degreeRadius) {
				    degreeRadius = newdegreeRadius;
                }
		    }
            gfile <<  " region "<< j << " degreeRadius  "<< degreeRadius << "  " <<endl << endl;
            regionCount++;
        } //end of if on non-zero regions
    } //end of loop on numpoints

//computing average values over all regions
    avDegreeAngle = 0;
	avDegreeRadius = 0;
	occupancy = 0;
    avAbsCorrelation = 0;
	tempPerimeter = 0;
	tempArea = 0;
	countRegions = 0;
	cout << "just before renewcorrelate_edges morph2 " << endl;
	M->random_correlate(max_comp, 2, Lzero);
	cout << "just after random correlate_edges morph2 " << endl;
    radiusDVector.resize(0);
    angleDVector.resize(0);
    angleVector.resize(0);
    radiusVector.resize(0);
    for (int j=0;j<numpoints;j++) {
	    if (M->regArea(j) != 0){
            countRegions++;
            avAbsCorrelation += M->renewcorrelate_edges(j,2,Lzero);
            occupancy += M->regNNfrac(j);
            cout << "after renNNfrac " << endl;
            tempArea = M->regArea(j);
            cout << " after regArea " << endl;
            tempPerimeter = M->regPerimeter(j);
            cout << " after regPerimeter " << endl;
            cout << "before sectorixe Dangle morph2 " << endl;
            angleDVector = M->sectorize_reg_Dangle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            cout << "after sectorixe Dangle morph2 " << endl;
            degreeAngle = L.find_zeroDAngle(angleDVector);
            avDegreeAngle += degreeAngle;
            cout << "after sectorixe Dangle morph2 " << endl;
            //radial degree
            degreeRadius = 0;
            radiusDVector = M->sectorize_reg_Dradius(j,numSectors, angleOffset, angleOffset + numSectors/2, S[j].NN);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            avDegreeRadius += degreeRadius;
            degfile3 << degreeAngle << " " << degreeRadius << " " << M->regNNfrac(j) << " " << tempArea << " " << tempPerimeter<<endl<<flush;
         } //end of if on non-zero regions
    } //end of loop on numpoints
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
	occupancy = occupancy / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
	jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation  <<endl;
    return 0;
};
