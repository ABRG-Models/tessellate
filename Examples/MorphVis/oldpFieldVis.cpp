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
#include <morph/Config.h>
#include <morph/Scale.h>
#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>

// Helper function to save PNG images with a suitable name
void savePngs (const std::string& logpath, const std::string& name,
               unsigned int frameN, morph::Visual& v)
{
    std::stringstream ff1;
    ff1 << logpath << "/" << name<< "_";
    ff1 << std::setw(5) << std::setfill('0') << frameN;
    ff1 << ".png";
    v.saveImage (ff1.str());
}

using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::cout;
using std::endl;
using std::runtime_error;
using namespace morph;
using namespace std;

#ifdef COMPILE_PLOTTING
class vtxVisual : public morph::VisualModel
{
public:
    // class data
    std::vector<std::vector<morph::Vector<FLT,3>>> vtx;
    unsigned int size;
    FLT linewidth = 0.00390625;
    FLT radiusFixed = 0.00390625;
    VBOint idx = 0;
    morph::Vector<FLT,3> uz = {1.0f, 1.0f, 1.0f};
    // class constructors
    vtxVisual(GLuint sp, GLuint tsp, std::vector<std::vector<morph::Vector<FLT, 3>>> & _vtx, FLT line, FLT radius)
    {
        this->shaderprog = sp;
        this->tshaderprog = tsp;
        this->vtx = _vtx;
        this->size = _vtx.size();
        this->linewidth = line;
        this->radiusFixed = radius;
    }
    //class methods
    void initializeVertices() {
        this->initv();
        cout << " in initializeVertices after call to initv " << endl;
    }

    // Draw vertices for the net's actual locations
    void initv() {
    // Discs at the net vertices
        cout << "entering intitv " << endl;
        morph::Vector<FLT,3> puckthick = { 0, 0, 0.002 };
        morph::Vector<FLT,3> linethick = {0, 0, 0.001};
        morph::Vector<FLT,3> colStart = {0.0, 0.0, 0.0};
        morph::Vector<FLT,3> colEnd = {0.0, 0.0, 0.0};
        for (unsigned int j=0; j < this->size; ++j) {
            unsigned int jsize = vtx[j].size();
            for (unsigned int i=0; i < jsize; i++) {
                cout << " initv vtx " << " i " << i << " j " << j << " is " << vtx[j][i] << endl;
                this->computeTube (this->idx,
                      this->vtx[j][i]+puckthick,
                      this->vtx[j][i]-puckthick,
                      morph::Vector<FLT,3>({1,0,0}), morph::Vector<FLT,3>({0,1,0}),
                      colStart, colEnd,
                      radiusFixed, 16);
            }
            cout << "after inner loop for computeTube " << endl;
        }
        cout << "after connect tube" << endl;
        // Connection lines
        for (unsigned int j=0; j < size; j++) {
            unsigned int jsize = vtx[j].size();
            for (unsigned int i=0; i < jsize; i++) {
                morph::Vector<FLT, 3> c1 = this->vtx[j][i] + linethick;
                morph::Vector<FLT, 3> c2 = this->vtx[j][(i+1)%jsize] + linethick;
                this->computeLine (idx, c1, c2, this->uz, colStart, colEnd, linewidth, linewidth/4);
            }
        }
    }
}; //end of class vtxVisual
#endif

int main (int argc, char **argv)
{
    if (argc < 2) {
      std::cout << " supply the json file as arg[1]" << endl;
      return -1;
    }
    string jsonfile = argv[1];
    string logpath = argv[2];
    //  open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
#ifdef SINGLE
        float dt = conf.getFloat("dt",0.001);
        float Dn = conf.getFloat("Dn",1.0);
        float Dchi = conf.getFloat("Dchi",0.0);
        float Dc = conf.getFloat("Dc",0.3);
        float xspan = conf.getFloat("xspan",5.0);
        float boundaryFalloffDist = conf.getFloat("boundaryFalloffDist",0.0078);
        float aNoiseGain = conf.getFloat("aNoiseGain",0.1);
        float nnInitialOffset = conf.getFloat("nnInitialOffet", 1.0);
        float ccInitialOffset = conf.getFloat("ccInitialOffset",2.5);
        float diffTol = conf.getFloat("diffTol",1e-8);
        float lengthScale = conf.getFloat("lengthScale",29.0f);
        float exponent = conf.getFloat("exponent", -100.0f);
        float radExp = conf.getFloat("radExp", -100.0f);
        float radMix = conf.getFloat("radMix", 1.0f);
#else
        double dt = conf.getDouble("dt",0.001);
        double Dn = conf.getDouble("Dn",1.0);
        double Dchi = conf.getDouble("Dchi",0.0);
        double Dc = conf.getDouble("Dc",0.3);
        double xspan = conf.getDouble("xspan",5.0);
        double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
        double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
        double nnInitialOffset = conf.getDouble("nnInitialOffet", 1.0);
        double ccInitialOffset = conf.getDouble("ccInitialOffset",2.5);
        double diffTol = conf.getDouble("diffTol",1e-8);
        double lengthScale = conf.getDouble("lengthScale",29.0);
        double exponent = conf.getDouble("exponent", -100.0);
        double radExp = conf.getDouble("radEXp", -100.0);
        double radMix = conf.getDouble("radMix", 1.0);
#endif
    int numSectors = conf.getInt("numSectors",24);
    int scale = conf.getInt("scale",8);
    unsigned int numsteps = conf.getUInt("numsteps",1000000);
    unsigned int numAdjust = conf.getUInt("numsteps",10000000);
    int plotEvery = conf.getInt("plotEvery",1000);
    int checkEvery = conf.getInt("checkEvery",1000);
    int fov = conf.getInt("fov",50);
    string iter = conf.getString("iter","0");
    bool LfixedSeed = conf.getBool("LfixedSeed",0);
    bool LDn = conf.getBool("LDn",false);
    //bool overwrite_logs = conf.getBool("overwrite_logs",true);
    bool skipMorph  = conf.getBool("skipMorph",false);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    bool lHomogen = conf.getBool("lHomogen",false);
    int  iPolygon = conf.getInt("iPolygon", 2);
    int  iBoundZero = conf.getInt("iBoundZero", 2);
    unsigned int numpoints = conf.getInt("numpoints",41);
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << " iBoundZero " << iBoundZero << std::endl;
    cout << "logpath " << logpath << std::endl;
    ofstream afile (logpath + "/centroids.out",ios::app);
#ifdef COMPILE_PLOTTING
    vtxVisual* cv;
#endif


    // A ra2yyndo2yym uniform generator returning real/FLTing point types
    ofstream gfile ( logpath +  "/edges.out");
    ofstream jfile ( logpath + "/results.txt",ios::app);
    ofstream degfile1 (logpath + "/degree1.data", ios::app);
    ofstream degfile2 (logpath + "/degree2.data", ios::app);
    ofstream degfile3 (logpath + "/degree3.data", ios::app);
    if (!gfile) {
        std::cout << "error opening gfile" << std::endl;
    }
    if (!jfile) {
        std::cout << "error opening jfile" << std::endl;
    }
    if (!degfile1) {
        std::cout << "error opening degfile1" << std::endl;
    }
    if (!degfile2) {
        std::cout << "error opening degfile2" << std::endl;
    }
    if (!degfile3) {
        std::cout << "error opening degfile3" << std::endl;
    }
    //scale time step according to Dn
    dt = dt / Dn;
    //set seed
    std::cout << "LfixedSeed = " << LfixedSeed << std::endl;
    unsigned int seed;
    chrono::milliseconds ms1 = duration_cast<milliseconds>(system_clock::now().time_since_epoch());
    if (LfixedSeed)
        seed = 596927;
    else
        seed = static_cast<unsigned int> (ms1.count());

    morph::RandUniform<FLT> ruf(seed);
// initialise DRegion class setting scale
    DRegion M(scale,xspan,logpath,numpoints,iPolygon); //create tessellation
    M.setCreg(); //set counts to identify inner boundaries
    cout << "after setCreg" << std::endl;
    M.setInternalBoundary(); //set internal boundaries
    cout << "after setInternalBoundary" << std::endl;
    //next line only needed for rectangular type domains
    if (iPolygon == 0) {
        M.cornerVertices(); //sets the rectangle corner hexes as vertices
    }
    cout << "before dissect_boundary " << endl;
    vector<std::pair<FLT,FLT>> cGravity;
    cGravity = M.dissectBoundary(); //dissect region boundary
    cout << "Edges size = " << M.edges.size() << endl;
	M.setRadialSegments(); //set the radial segments for regions
    int inReg = 0;
    inReg = M.setInnerRegion(); //set mask array for inner regions
    int rC = 0;
    //check for correct number of inner regions
    for (unsigned int j=0; j<numpoints;j++) {
        if (!M.innerRegion[j]) rC++;
    }
    if (inReg != rC) {
        cout << "Error: setInnerRegion returns " << inReg << " but outerRegions " << rC << endl;
        //return -1;
    }
    else {
        cout << "Success: setInnerRegion no of outer regions " << inReg << endl;
    }
    cout << "after first setRadialSegments " << endl;
// include the analysis methods
    Analysis L;

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate ksSolver
     */
    vector<ksSolver> S;
    /*
     * now we set the arrays for variable Dn, Dchi and Dc
     */
    vector<FLT> DnVal;
    DnVal.resize(numpoints,0.0);
    vector<FLT> DchiVal;
    DchiVal.resize(numpoints,0.0);
    vector<FLT> DcVal;
    DcVal.resize(numpoints,0.0);
    vector<FLT> morph0Area;
    morph0Area.resize(numpoints,0.0);
    //only need this for when we investigate random jitter in parameters
    /*
    for (unsigned int j = 0; j<numpoints; j++) {
        if (LDn) {
            DchiVal[j] = Dchi * (1.0 + ruf.get()*0.5) ;
            DnVal[j] = Dn *  (1.0 + ruf.get()*0.5);
            DcVal[j] = Dc * (1.0 + ruf.get()*0.5);
            cout << "j " << j << " DnVal[j] " <<  DnVal[j] <<  " DchiVal[j] " <<  DchiVal[j] <<  " DcVal[j] " << DcVal[j] << endl;
        }
        else {
            DchiVal[j] = Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }
    }
    */
    /*
     * Set boundaries based on polygons derived from the Dregion tesselation.
     * stored in curvedBoundary
     */
    M.populateBoundPolygon(1);
    cout << "just after setting polygonal boundaries " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j],lengthScale));
        cout << "in the loop populating the ksVector morph0 "<< j <<endl;
    }
    cout << "first setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
    // repopulate the regions
    cout << "just before populating the inner regions Morph 0" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewRegion(j,S[j].Hgrid->hexen);
    }
    cout << "just before populating the inner boundary Morph 0" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewBoundary(j,S[j].Hgrid->hexen);
    }
    cout << "just before calculating regionSize " << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << " Region size " << M.regionHex[j].size() << endl;
    }
	cout << "just after populating the regions from the ksSolver vector"<<endl;
    //clear global edges map
    M.edges_clear();
    cout << "just before renewDissect first time" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewDissect(j,0);
    }
    cout << "Edges size " << M.edges.size() << endl;
 #ifdef COMPILE_PLOTTING
// now draw the intial tesselation
    FLT hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;
    // Parameters from the config that apply only to plotting:
    // Should the plots be saved as png images?
    //const bool saveplots = conf.getBool ("saveplots", true);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    const bool vidframes = conf.getBool ("vidframes", true);
    unsigned int framecount = 0;

    // Window width and height
    const unsigned int win_width = conf.getUInt ("win_width", 2050UL);
    unsigned int win_height_default = static_cast<unsigned int>(0.88f * (float)win_width);
    const unsigned int win_height = conf.getUInt ("win_height", win_height_default);
    cout << "just before new Visual object" << endl;
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    morph::Visual * v1;
    v1 = new morph::Visual(win_width, win_height, "Tessellation0 ");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v1->backgroundWhite();
    //orthographic projection
    //v1->ptype = morph::perspective_type::orthographic;
    //v1->ortho_bl = {-1.0f, -1.0f};
    //v1->ortho_tr = {1.0f, 1.0f};
    // You can tweak the near and far clipping planes
    v1->zNear = 0.001;
    v1->zFar = 100;
    // And the field of view of the visual scene.
    v1->fov = fov;
    // You can lock movement of the scene
    v1->sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v1->setZDefault (conf.getFloat ("z_default", 0.0f));
    v1->setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v1->scenetrans_stepsize = 0.5;
    cout << "end of setting new Visual object" << endl;
    //make it current
    v1->setCurrent();

    // if using plotting, then set up the render clock
    steady_clock::time_point lastrender = steady_clock::now();
/*
// to draw the tessellation
// first convert M.vCoords into a morph::Vect
    cout << "just beofere creating vtxVector" << endl;
    std::vector<std::vector<morph::Vector<FLT,3>>> vtxVector;
    vtxVector.resize(M.vCoords.size());
    cout << "just after creating vtxVector" << endl;
    for (unsigned int j=0;j<numpoints;j++) {
        for (unsigned idx = 0; idx<M.vCoords[j].size(); idx++) {
            vtxVector[j].push_back(M.hGeo->point2vect3(M.vCoords[j][idx]));
        }
    }
    cout << "after filling vtxVector " << endl;
    // now instantiate vtxVisual
    cv = new vtxVisual(v1->shaderprog, v1->tshaderprog, vtxVector, M.ds, M.ds);
    cout << " after creating vtxVisual " << endl;
    cv->finalize();
    cout << " after vtxVisual finalize " << endl;
    cv->addLabel("Tessellation morph 0", {0.4f, 1.5, 0.0f});
    v1->addVisualModel(cv);
    cout << "before rendering v1 " << endl;
    v1->render();
    cout << "after rendering v1 " << endl;
    std::stringstream frame;
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    frame << framenum++;
    frame << ".png";
    cout << " before save image " << endl;
    savePngs (logpath, "tessellation0", 0, *v1);
*/

#endif

// initialise the fields
    string fname = logpath + "/first.h5";
    cout << "just before first data read morph 0 "<< " lcontinue " << Lcontinue << " exponent " << exponent << endl;
    // set up the boundaryFade vector and compute distance to boundary
    // Note the calls must be in this order, boundaryDist set first
    for (unsigned int j=0;j<numpoints;j++) {
        S[j].Hgrid->computeDistanceToBoundary();
        S[j].setBoundaryFade(exponent, boundaryFalloffDist);
    }
// initialise with random field
    if (Lcontinue) {
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++) {
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
    else {
        for (unsigned int j=0;j<numpoints;j++) {
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                FLT choice = ruf.get();
                if (choice > 0.5) {
                    if (lHomogen == true) {
                        S[j].NN[h.vi] = - ruf.get() * aNoiseGain   + nnInitialOffset;
                        S[j].CC[h.vi] = - ruf.get() * aNoiseGain   + ccInitialOffset;
                    }
                    else {
                        S[j].NN[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                        S[j].CC[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                    }
                }
                else {
                    if (lHomogen == true) {
                        S[j].NN[h.vi] = ruf.get() * aNoiseGain  + nnInitialOffset;
                        S[j].CC[h.vi] = ruf.get() * aNoiseGain  + ccInitialOffset;
                    }
                    else {
                        S[j].NN[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                        S[j].CC[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                    }
                }
                //what about the boundary?
                if (iBoundZero == 0) {
                    if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                        S[j].NN[h.vi] = S[j].NN[h.vi] * S[j].boundaryFade[h.vi];
                        S[j].CC[h.vi] = S[j].CC[h.vi] * S[j].boundaryFade[h.vi];
                    }
                }
                else if (iBoundZero == 1) {
                    if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                        S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * S[j].boundaryFade[h.vi] + nnInitialOffset;
                        S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * S[j].boundaryFade[h.vi] + ccInitialOffset;
                    }
                }
                else {
                    ;
                }
            } //end of loop over single region
        }//end of loop over regions
    } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

 #ifdef COMPILE_PLOTTING
    // Spatial offset, for positioning of HexGridVisuals
    Vector<float> spatOff;
    float xzero = 0.0f;

    // A. Offset in x direction to the left.
    xzero = 0.5 * M.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    Scale<FLT,float> zscale; zscale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    Scale<FLT,float> cscale; cscale.do_autoscale = true;

    unsigned int Agrid[numpoints];

    v1->addLabel("Laplacian", {0.4f, 1.3f, 0.0f});
    v1->addLabel("Field", {-0.4f, 1.3f, 0.0f});
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Agrid[j] = v1->addVisualModel (new HexGridVisual<FLT> (v1->shaderprog,
                                                                    v1->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::GreyscaleInv));
#ifdef RANDOM
        }//end of loop on inner regions
#endif
    }//end of loop over regions

    // A. Offset in x direction to the left.
    xzero -= 1.0 * M.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    unsigned int Bgrid[numpoints];
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Bgrid[j] = v1->addVisualModel (new HexGridVisual<FLT> (v1->shaderprog,
                                                                    v1->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::GreyscaleInv));
#ifdef RANDOM
        }//end of loop on inner regions
#endif
    }//end of loop over regions
#endif
    // begin time stepping loop unmorphed grid solved via schSolver
    // set up vectors for determining convergence
    std::vector<FLT> NNdiff;
    std::vector<std::vector<FLT>> NNpre;
    std::vector<std::vector<FLT>> NNcurr;
    NNdiff.resize(numpoints);
    NNpre.resize(numpoints);
    NNcurr.resize(numpoints);
    FLT NNdiffSum = 0.0f;
  //initilise all Apre vectors to above possible field
    for (unsigned int j=0; j<numpoints;j++) {
        NNpre[j].resize(S[j].NN.size(),10000.0);
    }
    // begin morph0 time stepping loop
    for (unsigned int i=0;i<numsteps;i++) {
   	for (unsigned int j = 0;j<numpoints;j++) { //loop over all regions
            S[j].stepEuler(dt, Dn, Dchi, Dc);
            if (i%checkEvery == 0) {
                NNdiffSum = 0.0;
                NNcurr[j] = S[j].NN;
                NNdiff[j] = L.normedDiff(NNpre[j], NNcurr[j]);
                NNpre[j] = NNcurr[j];
                NNdiffSum += fabs(NNdiff[j]);
            }
            /*
            if (i%numAdjust == 0) {
                S[j].signal(radExp, radMix, aNoiseGain);
            }
            */
        } //end of loop over regions
        if (i%checkEvery == 0) {
            cout << "NNdiffSum " << NNdiffSum << " i = " << i << endl;
            if (NNdiffSum/(numpoints*1.0) < diffTol) {
                cout << "morphed converged step " << i << " field diff " << NNdiffSum/(1.0*numpoints) << " diffTol " << diffTol << std::endl;
                break;
            }
        }
#ifdef COMPILE_PLOTTING
        if ((i % plotEvery) == 0) {
            cout << "step " << i << " reached" << endl;
            for (unsigned int j=0; j<numpoints;j++) {
#ifdef RANDOM
                if (M.innerRegion[j]) { //only display inner regions
#endif
                    vector<FLT> regionNN;
                    regionNN = L.normalise(S[j].lapNN);
                    VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v1->getVisualModel (Agrid[j]);
                    avm->updateData (&regionNN);
                    avm->clearAutoscaleColour();
                    regionNN = L.normalise(S[j].NN);
                    VisualDataModel<FLT>* avm1 = (VisualDataModel<FLT>*)v1->getVisualModel (Bgrid[j]);
                    avm1->updateData (&regionNN);
                    avm1->clearAutoscaleColour();
                    v1->render();
                    if (i%plotEvery == 0) {
                        std::cout << std::setprecision(16) << " sum of NN " << S[j].sum_NN << " sum of lapNN " << S[j].sum_lapNN << " step " << i << " plotEvery " << plotEvery << std::endl;
                    }
#ifdef RANDOM
                }//end of if on inner regions
#endif
            }
            if (vidframes) {
                savePngs (logpath, "nn0", framecount, *v1);
                ++framecount;
            }
            else {
                savePngs (logpath, "nn0", i , *v1);
            }
        }
#endif
    }//end of time stepping
    //code run at end of timestepping
    //first save the  ofstream outFile;
    cout << "just before first data write morph 0" << endl;
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
    }
    fdata.add_val ("/Dchi", Dchi);
    fdata.add_val ("/Dn", Dn);
    fdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;
    gfile << endl << "analysis on first morphing iteration " << endl;

    //declaration of variables needed for analysis
    vector <FLT> angleVector;
    vector <FLT> radiusVector;
    vector <int> radiusDVector;
    int degreeRadius;
    int degreeAngle;
    FLT tempArea = 0.0;
    FLT tempPerimeter = 0.0;
    int angleOffset = 0;
    int radiusOffset = 0;
    FLT avDegreeRadius = 0.0;
    FLT avDegreeAngle = 0;
    FLT occupancy = 0.0;
    int countRegions = 0;
    FLT avAbsCorrelation = 0.0;
    const int max_comp = numpoints*3;
    for (unsigned int j=0;j<numpoints;j++) {
        M.NN[j] = S[j].NN; //write the NN fields to the DRegion array
    }
    for (unsigned int j=0;j<numpoints;j++) {
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            countRegions++;
            occupancy += M.regNNfrac(j);
            tempArea = M.regArea(j);
            tempPerimeter = M.renewRegPerimeter(j);
            std::cout << "just before renewcorrelate edges morph0" << std::endl;
            avAbsCorrelation += M.renewcorrelate_edges(j,1);
            std::cout << "just before sectorize reg" << std::endl;
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            std::cout << "just before extrema angle" << std::endl;
            degreeAngle = L.find_extrema_angle(angleVector);
            avDegreeAngle += degreeAngle;
            std::cout << "just before sectorize radius" << std::endl;
            //radial degree
            degreeRadius = 0;
            //radiusVector = M.sectorize_reg_radius(j,numSectors/2, angleOffset, angleOffset + numSectors, S[j].NN);
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors/2, angleOffset, angleOffset + numSectors, S[j].NN);
            std::cout << "just before extrema radius morph0" << std::endl;
            //degreeRadius = L.find_extrema_radius(radiusVector);
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            //degreeRadius += 1; //there must always be an extremum at the boundary
            avDegreeRadius += degreeRadius;
            degfile1 << degreeAngle/2 << " " << degreeRadius  << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
#ifdef RANDOM
        } //end of if on non-zero regions
#endif
    } //end of loop on NUMPOINTs
    cout << "just before renewcorrelate_edges morph0 " << endl;
    //avAbsCorrelation = M.correlate_edges(0);
    M.random_correlate(max_comp, 1);
    cout << "just after randomcorrelate_edges morph0 " << endl;
    if (countRegions == 0) {
        cout << "Error zero regionss counted in second analysis morph 0" << endl;
        return -1;
    }
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
    occupancy = occupancy / (1.0 * countRegions);
    cout << "Regions counted in first analysis loop " << countRegions << endl;
    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;

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

    /*
     * now create an array of morphed regions, first morph
     */

    if (skipMorph) return 0;
    cout << "just before setting curved boundaries first morph" <<endl;
// section for solving on the curved boundaries
    S.resize(0);
    //set radius for creating circular regions also calculate the adjusted Dn values
	for (unsigned int j = 0;j<numpoints;j++) {
        morph0Area[j] = M.hexArea*M.regArea(j);
    }
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(1);
    cout << "just after setting curved boundaries morph 1 " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j],lengthScale));
    }
    for (unsigned int j = 0;j<numpoints;j++) {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j],lengthScale));
        cout << "in the loop populating the ksVector morph1   "<< j << " size S[j] " << S[j].Hgrid->hexen.size() << endl;
    }
    cout << "first morph  setting of centroids" << endl;
    cout << "first morph  setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
    // repopulate the regions
    cout << "just before populating the inner regions morph 1" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewRegion(j,S[j].Hgrid->hexen);
    }
    cout << "just before populating the inner boundary morph 1" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewBoundary(j,S[j].Hgrid->hexen);
    }
    cout << "second setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        afile << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
    }
	cout << "just after populating the ksVector"<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        FLT area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] = Dchi * (area / morph0Area[j]);
            DnVal[j] = Dn *  (area / morph0Area[j]);
            DcVal[j] = Dc * (area / morph0Area[j]);
        }
        else {
            DchiVal[j] = Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }

        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
// initialise the fields
    string gname = logpath + "/second.h5";
    cout<< "just before second data read"<< " lcontinue " << Lcontinue <<endl;
    // Note the calls must be in this order, boundaryDist set first
    for (unsigned int j=0;j<numpoints;j++) {
        S[j].Hgrid->computeDistanceToBoundary();
        S[j].setBoundaryFade(exponent, boundaryFalloffDist);
    }
// initialise with random field
    if (Lcontinue) {
        morph::HdfData ginput(gname,1);
        cout << "just after trying to open ../logs/second.h5" << endl;
        for (unsigned int j=0;j<numpoints;j++) {
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
    else {
        for (unsigned int j=0;j<numpoints;j++) {
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                FLT choice = ruf.get();
                if (choice > 0.5) {
                    if (lHomogen == true) {
                        S[j].NN[h.vi] = - ruf.get() * aNoiseGain   + nnInitialOffset;
                        S[j].CC[h.vi] = - ruf.get() * aNoiseGain   + ccInitialOffset;
                    }
                    else {
                        S[j].NN[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                        S[j].CC[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                    }
                }
                else {
                    if (lHomogen == true) {
                        S[j].NN[h.vi] = ruf.get() * aNoiseGain  + nnInitialOffset;
                        S[j].CC[h.vi] = ruf.get() * aNoiseGain  + ccInitialOffset;
                    }
                    else {
                        S[j].NN[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                        S[j].CC[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                    }
                }
                //what about the boundary?
                if (iBoundZero == 0) {
                    if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                        S[j].NN[h.vi] = S[j].NN[h.vi] * S[j].boundaryFade[h.vi];
                        S[j].CC[h.vi] = S[j].CC[h.vi] * S[j].boundaryFade[h.vi];
                    }
                }
                else if (iBoundZero == 1) {
                    if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                        S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * S[j].boundaryFade[h.vi] + nnInitialOffset;
                        S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * S[j].boundaryFade[h.vi] + ccInitialOffset;
                    }
                }
                else {
                    ;
                }
            } //end of loop over a single region
        }//end of loop over regions
    } //end of else on Lcontinue
    cout <<  "just after field creation first morph" << endl;
#ifdef COMPILE_PLOTTING
    morph::Visual * v2;
    //morph::Vector<FLT,3> offset ={1.0,0.0,0.0};
    v2 = new morph::Visual(win_width, win_height, "Tessellation1 ");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v2->backgroundWhite();
/*
    v2->ptype = morph::perspective_type::orthographic;
    v2->ortho_bl = {-1.0f, -1.0f};
    v2->ortho_tr = {1.0f, 1.0f};
*/
    // You can tweak the near and far clipping planes
    v2->zNear = 0.001;
    v2->zFar = 500;
    // And the field of view of the visual scene.
    v2->fov = fov;
    // You can lock movement of the scene
    v2->sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v2->setZDefault (conf.getFloat ("z_default", -5.0f));
    v2->setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v2->scenetrans_stepsize = 0.5;
    /*
    // now instantiate vtxVisual
    cv = new vtxVisual(v2->shaderprog, v2->tshaderprog, vtxVector, M.ds, M.ds);
    cout << " after creating vtxVisual " << endl;
    cv->finalize();
    cout << " after vtxVisual finalize " << endl;
    cv->addLabel("Tessellation morph 1", {0.4f, 1.4f, 0.0f});
    v2->addVisualModel(cv);
    cout << "before rendering v2 " << endl;
    v2->render();
    cout << "after rendering v2 " << endl;
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    frame << framenum++;
    frame << ".png";
    cout << " before save image " << endl;
    savePngs (logpath, "tessellation1", 0, *v2);
    v2->setCurrent();
*/
    // A. Offset in x direction to the left.
    xzero = 0.5 * M.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    //Scale<FLT,float> zscale; zscale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    //Scale<FLT,float> cscale; cscale.do_autoscale = true;

    unsigned int Cgrid[numpoints];

    v2->addLabel("Laplacian", {0.4f, 1.3f, 0.0f});
    v2->addLabel("Field", {-0.4f, 1.3f, 0.0f});
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Cgrid[j] = v2->addVisualModel (new HexGridVisual<FLT> (v2->shaderprog,
                                                                    v2->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::GreyscaleInv));
#ifdef RANDOM
        }//end of loop on inner regions
#endif
    }//end of loop over regions

    // A. Offset in x direction to the left.
    xzero -= 1.0 * M.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    unsigned int Dgrid[numpoints];
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Dgrid[j] = v2->addVisualModel (new HexGridVisual<FLT> (v2->shaderprog,
                                                                    v2->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::GreyscaleInv));
#ifdef RANDOM
        }//end of loop on inner regions
#endif
    }//end of loop over regions
#endif
    // begin second time stepping loop after first morph
    std::cout << "before time stepping loop morph 1" << std::endl;
    NNdiff.resize(numpoints);
    NNpre.resize(numpoints);
    NNcurr.resize(numpoints);
    NNdiffSum = 0.0;

    //initilise all Apre vectors to above possible field
    for (unsigned int j=0; j<numpoints;j++) {
        NNpre[j].resize(S[j].NN.size(),1000.0);
    }
#ifdef COMPILE_PLOTTING
    //reset framecount
    framecount = 0;
#endif
    //start of time-stepping loo
    for (unsigned int i=0;i<numsteps;i++) {
        for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
            S[j].step(dt, Dn, Dchi, Dc);
            if (i%checkEvery == 0) {
                NNdiffSum = 0.0;
                NNcurr[j] = S[j].NN;
                NNdiff[j] = L.normedDiff(NNpre[j], NNcurr[j]);
                NNpre[j] = NNcurr[j];
                NNdiffSum += NNdiff[j];
            }
            /*
            if (i%numAdjust == 0 ) {
                S[j].signal(radExp, radMix, aNoiseGain);
            }
            */
        }
        if (i%checkEvery == 0) {
            cout << "NNdiffSum " << NNdiffSum << " i = " << i << endl;
            if (NNdiffSum/(numpoints*1.0) < diffTol) {
                cout << "morphed converged step " << i << " field diff " << NNdiffSum/(1.0*numpoints) << " diffTol " << diffTol << std::endl;
                break;
            }
        }

#ifdef COMPILE_PLOTTING
        if ((i % plotEvery) == 0) {
            std::cout << "in plotEvery of time loop i " << i << std::endl;
            for (unsigned int j=0; j<numpoints;j++) {
#ifdef RANDOM
                if (M.innerRegion[j]) { //only display inner regions
#endif
                    vector<FLT> regionNN;
                    regionNN = L.normalise(S[j].lapNN);
                    VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v2->getVisualModel (Cgrid[j]);
                    avm->updateData (&regionNN);
                    avm->clearAutoscaleColour();
                    regionNN = L.normalise(S[j].NN);
                    VisualDataModel<FLT>* avm1 = (VisualDataModel<FLT>*)v2->getVisualModel (Dgrid[j]);
                    avm1->updateData (&regionNN);
                    avm1->clearAutoscaleColour();
                    if (i%plotEvery == 0) {
                        std::cout << std::setprecision(16) << " sum of NN " << S[j].sum_NN << " sum of lapNN " << S[j].sum_lapNN << " step " << i << " plotEvery " << plotEvery << std::endl;
                    }
    #ifdef RANDOM
                }
    #endif
            }
            v2->render();
            if (vidframes) {
                savePngs (logpath, "nn1", framecount, *v2);
                ++framecount;
            }
            else {
                savePngs (logpath, "nn1", i , *v2);
            }
        }
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v2.render().
        //steady_clock::duration sincerender = steady_clock::now() - lastrender;
        //if (duration_cast<milliseconds>(sincerender).count() > 17000) { // 17 is about 60 Hz
        //    glfwPollEvents();
        //    v2->render();
        //    lastrender = steady_clock::now();
        //}
#endif
    } //end of numsteps loop
    //code run at end of timestepping
    //first save the  ofstream outFile;
    cout << "just before second data read " << endl;
    morph::HdfData gdata(gname);
    for (unsigned int j=0;j<numpoints;j++) {
        std::string nstr = "n" + to_string(j);
        char * nst = new char[nstr.length()+1];
        std::strcpy(nst,nstr.c_str());
        std::string ccstr = "c" + to_string(j);
        char * ccst = new char[ccstr.length()+1];
        std::strcpy(ccst,ccstr.c_str());
        cout << "labels "<< nst <<" , " << nstr <<","<< ccst<< "," << ccstr <<endl;
        gdata.add_contained_vals(ccst,S[j].CC);
        gdata.add_contained_vals(nst,S[j].NN);
     }
     gdata.add_val ("/Dchi", Dchi);
     gdata.add_val ("/Dn", Dn);
     gdata.add_val ("/Dc",Dc);
     cout << " just after writing data "  << endl;
     // write the NN and CC vals for each region
      gfile << endl << "analysis on first morphing iteration " << endl;
    angleVector.resize(0);
    radiusDVector.resize(0);
    radiusVector.resize(0);


    M.edges_clear();
    // swap the radialAngles to the mCoords
    // M.swapRadialSegments(false);
    // redissect the boundaries
    cout << "just before renewDissect second time" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewDissect(j,1);
    }
    cout << "Edges size " << M.edges.size() << endl;


    avDegreeAngle = 0;
    avDegreeRadius = 0;
    occupancy = 0;
    countRegions = 0;
    tempArea = 0;
    tempPerimeter = 0;
    avAbsCorrelation = 0;
    M.random_correlate(max_comp,2);
    cout << "just after randomcorrelate_edges morph1 " << endl;
    for (unsigned int j=0;j<numpoints;j++) {
        M.NN[j] = S[j].NN; //write the NN fields to the DRegion array
    }
    for (unsigned int j=0;j<numpoints;j++) {
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            countRegions++;
            std::cout << "just before occupancy morph1" << std::endl;
            occupancy += M.regNNfrac(j);
            std::cout << "just before regArea morph1" << std::endl;
            tempArea = M.regArea(j);
            std::cout << "just before renewPerimeter morph1" << std::endl;
            tempPerimeter = M.renewRegPerimeter(j);
            std::cout << "just before renewcorrelate_edges morph1" << std::endl;
            avAbsCorrelation += M.renewcorrelate_edges(j,2);
            std::cout << "just before sectorize_reg_angle morph1" << std::endl;
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            std::cout << "just before find extrema angle morph1" << std::endl;
            degreeAngle = L.find_extrema_angle(angleVector);
            avDegreeAngle += degreeAngle;
            //radial degree
            degreeRadius = 0;
            std::cout << "just before sectorize_reg_Dradius morph1" << std::endl;
            //radiusVector = M.sectorize_reg_radius(j,numSectors/2, angleOffset, angleOffset + numSectors, S[j].NN);
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors/2, angleOffset, angleOffset + numSectors, S[j].NN);
            std::cout << "just before find extrema radius morph1" << std::endl;
            degreeRadius = L.find_zeroDRadius(radiusDVector);
            //degreeRadius = L.find_extrema_radius(radiusVector);
            //degreeRadius += 1; //there is always an extremum at the boundary
            avDegreeRadius += degreeRadius;
            degfile2 << degreeAngle/2 << " " << degreeRadius << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
    #ifdef RANDOM
        } //end of if on non-zero regions
    #endif
    } //end of loop on NUMPOINTs
    cout << "just before renewcorrelate_edges morph1 " << endl;
    //avAbsCorrelation = M.correlate_edges(0);
    M.random_correlate(max_comp, 1);
    cout << "just after randomcorrelate_edges morph1 " << endl;
    if (countRegions == 0) {
        cout << "Error zero regionss counted in third analysis morph 1" << endl;
        return -1;
    }
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation / (1.0 * countRegions);
    occupancy = occupancy / (1.0 * countRegions);
    jfile <<Dn<< " "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation << endl;
//end of integration after first morphing
//begin second morphing
    if (skipMorph) return 0;
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(0);
    S.resize(0);
    cout << "just after setting curved boundaries morph 2 " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++)
    {
        S.push_back(ksSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j], lengthScale));
        cout << "in the loop populating the ksVector morph2 "<< j <<endl;
    }
    cout << "just after populating the ksVector"<<endl;
    // repopulate the regions
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewRegion(j,S[j].Hgrid->hexen);
    }
    cout << "after repopulate regions morph2" << endl;
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewBoundary(j,S[j].Hgrid->hexen);
        //M.renewCentroids(j);
    }
    cout << "after repopulate boundary morph2" << endl;
//    swap the radialAngles to the mCoords
    M.edges_clear();
    M.swapRadialSegments(false);
    // redissect the boundaries
    for (unsigned int j=0;j<numpoints;j++)
    {
            M.renewDissect(j,2);
    }
    cout << "Edges size " << M.edges.size() << endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        FLT area = M.hexArea*M.regArea(j);
        if (LDn) {
            DchiVal[j] =  Dchi * (area / morph0Area[j]);
            DnVal[j] = Dn * (area / morph0Area[j]);
            DcVal[j] = Dc * (area / morph0Area[j]);
        }
        else {
            DchiVal[j] =  Dchi;
            DnVal[j] = Dn;
            DcVal[j] = Dc;
        }
        cout << "DnVal region " << j << " = " << DnVal[j] << " Dn = " << Dn << " PI " << PI << endl;
    }
// now draw the intial tesselation
// initialise the fields
        cout<< "just before third data read"<< " lcontinue " << Lcontinue <<endl;
    // Note the calls must be in this order, boundaryDist set first
    for (unsigned int j=0;j<numpoints;j++) {
        S[j].Hgrid->computeDistanceToBoundary();
        S[j].setBoundaryFade(exponent, boundaryFalloffDist);
    }
// initialise with random field
    string hname = logpath + "/third.h5";
    if (Lcontinue) {
        morph::HdfData hinput (hname,1);
        for (unsigned int j=0;j<numpoints;j++){
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
        for (unsigned int j=0;j<numpoints;j++) {
            for (auto h : S[j].Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
                FLT choice = ruf.get();
                if (choice > 0.5) {
                    if (lHomogen == true) {
                        S[j].NN[h.vi] = - ruf.get() * aNoiseGain   + nnInitialOffset;
                        S[j].CC[h.vi] = - ruf.get() * aNoiseGain   + ccInitialOffset;
                    }
                    else {
                        S[j].NN[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                        S[j].CC[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                    }
                }
                else {
                    if (lHomogen == true) {
                        S[j].NN[h.vi] = ruf.get() * aNoiseGain  + nnInitialOffset;
                        S[j].CC[h.vi] = ruf.get() * aNoiseGain  + ccInitialOffset;
                    }
                    else {
                        S[j].NN[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                        S[j].CC[h.vi] = aNoiseGain * (exp(radExp * h.r*h.r)*radMix + ruf.get());
                    }
                }
                //what about the boundary?
                if (iBoundZero == 0) {
                    if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                        S[j].NN[h.vi] = S[j].NN[h.vi] * S[j].boundaryFade[h.vi];
                        S[j].CC[h.vi] = S[j].CC[h.vi] * S[j].boundaryFade[h.vi];
                    }
                }
                else if (iBoundZero == 1) {
                    if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                        S[j].NN[h.vi] = (S[j].NN[h.vi] - nnInitialOffset) * S[j].boundaryFade[h.vi] + nnInitialOffset;
                        S[j].CC[h.vi] = (S[j].CC[h.vi] - ccInitialOffset) * S[j].boundaryFade[h.vi] + ccInitialOffset;
                    }
                }
                else {
                    ;
                }
            } //end of loop over single region
        }//end of loop over regions
    } //end of else on Lcontinue
#ifdef COMPILE_PLOTTING
    morph::Visual * v3;
    v3 = new morph::Visual(win_width, win_height, "Tessellation2");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v3->backgroundWhite();
    // You can tweak the near and far clipping planes
    v3->zNear = 0.001;
    v3->zFar = 500;
    // And the field of view of the visual scene.
    v3->fov = fov;
    // You can lock movement of the scene
    v3->sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v3->setZDefault (conf.getFloat ("z_default", 0.0f));
    v3->setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v3->scenetrans_stepsize = 0.5;
    cout << "end of setting new Visual object" << endl;
    /*
    // now instantiate vtxVisual
    cv = new vtxVisual(v3->shaderprog, v3->tshaderprog, vtxVector, M.ds, M.ds);
    cout << " after creating vtxVisual " << endl;
    cv->finalize();
    cout << " after vtxVisual finalize " << endl;
    cv->addLabel("Tessellation morph 2", {0.4f, 1.4f, 0.0f});
    v3->addVisualModel(cv);
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    frame << framenum++;
    frame << ".png";
    v3->setCurrent();
    v3->render();
*/
    // A. Offset in x direction to the left.
    xzero = 0.5 * M.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    //zscale; zscale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    //Scale<FLT,float> cscale; cscale.do_autoscale = true;

    unsigned int Egrid[numpoints];

    v3->addLabel("Laplacian", {0.4f, 1.3f, 0.0f});
    v3->addLabel("Field", {-0.4f, 1.3f, 0.0f});
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Egrid[j] = v3->addVisualModel (new HexGridVisual<FLT> (v3->shaderprog,
                                                                    v3->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::GreyscaleInv));
#ifdef RANDOM
        }//end of loop on inner regions
#endif
    }//end of loop over regions

    // A. Offset in x direction to the left.
    xzero -= 1.0 * M.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    unsigned int Fgrid[numpoints];
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
#ifdef RANDOM
        if (M.innerRegion[j]) {
#endif
            vector<FLT> regionNN;
            //normalise over the region t
            regionNN = L.normalise(S[j].NN);
            Fgrid[j] = v3->addVisualModel (new HexGridVisual<FLT> (v3->shaderprog,
                                                                    v3->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionNN),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::GreyscaleInv));
#ifdef RANDOM
        }//end of loop on inner regions
#endif
    }//end of loop over regions
#endif

     //begin fourth time stepping loop after second morph
    NNdiff.resize(numpoints);
    NNpre.resize(numpoints);
    NNcurr.resize(numpoints);
    NNdiffSum = 0.0;

    //initilise all Apre vectors to above possible field
    for (unsigned int j=0; j<numpoints;j++) {
        NNpre[j].resize(S[j].NN.size(),1000.0);
    }
#ifdef COMPILE_PLOTTING
    //reset framecount
    framecount = 0;
#endif

    for (unsigned int i=0;i<numsteps;i++) {
        for (unsigned int j = 0;j<numpoints;j++) {
            S[j].step(dt, Dn, Dchi, Dc);
            if (i%checkEvery == 0) {
                NNdiffSum = 0.0;
                NNcurr[j] = S[j].NN;
                NNdiff[j] = L.normedDiff(NNpre[j], NNcurr[j]);
                NNpre[j] = NNcurr[j];
                NNdiffSum += NNdiff[j];
            }
            /*
            if (i%numAdjust == 0) {
                S[j].signal(radExp, radMix, aNoiseGain);
            }
            */
        }
        if (i%checkEvery == 0) {
            cout << "NNdiffSum " << NNdiffSum << " i = " << i << endl;
            if (NNdiffSum/(numpoints*1.0) < diffTol) {
                cout << "morphed converged step " << i << " field diff " << NNdiffSum/(1.0*numpoints) << " diffTol " << diffTol << std::endl;
                break;
            }
        }
#ifdef COMPILE_PLOTTING
    if((i % plotEvery) == 0) {
        for (unsigned int j=0; j<numpoints;j++) {
#ifdef RANDOM
            if (M.innerRegion[j]) { //only display inner regions
#endif
                vector<FLT> regionNN;
                regionNN = L.normalise(S[j].lapNN);
                VisualDataModel<FLT>* avm = (VisualDataModel<FLT>*)v3->getVisualModel (Egrid[j]);
                avm->updateData (&regionNN);
                avm->clearAutoscaleColour();
                regionNN = L.normalise(S[j].NN);
                VisualDataModel<FLT>* avm1 = (VisualDataModel<FLT>*)v3->getVisualModel (Fgrid[j]);
                avm1->updateData (&regionNN);
                avm1->clearAutoscaleColour();
                if (i%plotEvery == 0) {
                    std::cout << std::setprecision(16) << " sum of NN " << S[j].sum_NN << " sum of lapNN " << S[j].sum_lapNN << " step " << i << " plotEvery " << plotEvery << std::endl;
                }
#ifdef RANDOM
            }
#endif
        }
   	// rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v3.render().
        cout << " before rendering the graphics " << endl;
        v3->render();
        if (vidframes) {
            savePngs (logpath, "nn2", framecount, *v3);
            ++framecount;
        }
        else {
            savePngs (logpath, "nn2", i , *v3);
        }
        cout << " after rendering the graphics " << endl;
    }

#endif
} //end of numsteps loop

    //
       //nndisp.closeDisplay();
    //code run at end of timestepping
    cout << "just before the fourth data write" << endl;
    //first save the  ofstream outFile;
    morph::HdfData hdata(hname);
    cout <<"after creating hdata "<< hname << endl;
    for (unsigned int j=0;j<numpoints;j++) {
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
    //data.add_contained_vals("X",M.X[0]);
    //data.add_contained_vals("Y",M.X[1]);
    // hdata.add_val ("/Dchi", Dchi);
    // hdata.add_val ("/Dn", Dn);
    // hdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;

    gfile << endl << "analysis on second morphing iteration " << endl;
    angleVector.resize(0);
    radiusVector.resize(0);
    radiusDVector.resize(0);

//computing average values over all regions
    avDegreeAngle = 0;
    avDegreeRadius = 0;
    occupancy = 0;
    avAbsCorrelation = 0;
    tempPerimeter = 0;
    tempArea = 0;
    countRegions = 0;
    M.random_correlate(max_comp, 3);
    cout << "just after random correlate_edges morph2 " << endl;
    angleVector.resize(0);
    radiusVector.resize(0);
    for (unsigned int j=0;j<numpoints;j++) {
        M.NN[j] = S[j].NN; //write the NN fields to the DRegion array
    }
    for (unsigned int j=0;j<numpoints;j++) {
#ifdef RANDOM
        if (M.innerRegion[j]){
#endif
            countRegions++;
            std::cout << "just before occupancy morph2" << std::endl;
            occupancy += M.regNNfrac(j);
            std::cout << "just before regArea morph2" << std::endl;
            tempArea = M.regArea(j);
            std::cout << "just before renewPerimeter morph2" << std::endl;
            tempPerimeter = M.renewRegPerimeter(j);
            std::cout << "just before renewcorrelate_edges morph2" << std::endl;
            avAbsCorrelation += M.renewcorrelate_edges(j,2);
            std::cout << "just before sectorize_reg_angle morph2" << std::endl;
            angleVector = M.sectorize_reg_angle(j,numSectors,radiusOffset, numSectors, S[j].NN);
            std::cout << "just before find extrema angle morph2" << std::endl;
            degreeAngle = L.find_extrema_angle(angleVector);
            avDegreeAngle += degreeAngle;
            //radial degree
            degreeRadius = 0;
            std::cout << "just before sectorize_reg_Dradius morph2" << std::endl;
            //radiusVector = M.sectorize_reg_radius(j,numSectors/2, angleOffset, angleOffset + numSectors, S[j].NN);
            radiusDVector = M.sectorize_reg_Dradius(j,numSectors/2, angleOffset, angleOffset + numSectors, S[j].NN);
            std::cout << "just before find_extrema_radius morph2" << std::endl;
            degreeRadius = L.find_zeroDRadius(radiusDVector);
           // degreeRadius = L.find_extrema_radius(radiusVector);
           // degreeRadius += 1; //there must always be an extremum at the boundary
            avDegreeRadius += degreeRadius;
            degfile3 << degreeAngle/2 << " " << degreeRadius  << " " << M.regNNfrac(j) << " " << tempArea << " "<< tempPerimeter<<endl<<flush;
#ifdef RANDOM
         } //end of if on non-zero regions
#endif
    } //end of loop on NUMPOINTs
    cout << "just before renewcorrelate_edges morph2 " << endl;
    //avAbsCorrelation = M.correlate_edges(0);
    M.random_correlate(max_comp, 1);
    cout << "just after randomcorrelate_edges morph2 " << endl;
    if (countRegions == 0) {
        cout << "Error zero regionss counted in third analysis morph 2" << endl;
        return -1;
    }
    avDegreeAngle = avDegreeAngle / (1.0 * countRegions);
    avDegreeRadius = avDegreeRadius / (1.0 * countRegions);
    occupancy = occupancy / (1.0 * countRegions);
    avAbsCorrelation = avAbsCorrelation/(1.0 * countRegions);
    jfile <<Dn<<"  "<<Dchi<<" "<<Dc<<" "<<avDegreeAngle<<" "<<avDegreeRadius<<" "<<occupancy<<" "<<avAbsCorrelation  <<endl;
    cout << "end of program reached successfully!" << endl;
    return 0;
};
