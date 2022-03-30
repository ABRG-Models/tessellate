/*!
 * If COMPILE_PLOTTING is defined at compile time, then include the display and
 * plotting code. I usually put all the plotting code inside #defines like this so
 * that I can compile a version of the binary without plotting, for parameter searches
 * in which I am only going to be saving out HDF5 data.
 */
#include "shRegion.h"
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <morph/ColourMap.h>
#include <morph/VisualDataModel.h>
#include <morph/Scale.h>
#include "shCSolver.h"
#include "analysis.h"
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
#include <chrono>
#include <complex>

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
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HexGrid;
using morph::HdfData;
using morph::Tools;
using namespace morph;
using namespace std;
using namespace std::chrono;

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
      std::cout << "not enough arguments" << argc << endl;
      return -1;
    }
    string jsonfile = argv[1];
    //open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
 #ifdef SINGLE
        float dt = conf.getFloat("dt",0.0001);
        float epsilon = conf.getFloat("epsilon",5.0);
        float  g = conf.getFloat("g",5.0);
        float  xspan = conf.getFloat("xspan",5.0);
        float  boundaryFalloffDist = conf.getFloat("boundaryFalloffDist",0.0078);
        float  aNoiseGain = conf.getFloat("aNoiseGain",0.1);
        float  nnInitialOffset = conf.getFloat("nnInitialOffset",1.0);
        float  ccInitialOffset = conf.getFloat("ccInitialOffset", 2.5);
        float  x_default = conf.getFloat("x_default",0.0);
        float  y_default = conf.getFloat("y_default",0.0);
        float  wratio = conf.getFloat("wratio", 0.844f);
#else
        double dt = conf.getDouble("dt",0.0001);
        double epsilon = conf.getDouble("epsilon",5.0);
        double g = conf.getDouble("g",5.0);
        double xspan = conf.getDouble("xspan",5.0);
        double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
        double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
        double nnInitialOffset = conf.getDouble("nnInitialOffset",1.0);
        double ccInitialOffset = conf.getDouble("ccInitialOffset", 2.5);
        double  x_default = conf.getDouble("x_default",0.0);
        double  y_default = conf.getDouble("y_default",0.0);
        double  wratio = conf.getDouble("wratio", 0.844f);
#endif
    int scale = conf.getInt("scale",8);
    int numsteps = conf.getInt("numsteps",100);
    int numAdjust = conf.getInt("numAdjust",1000000);
    int numprint = conf.getInt("numprint",95);
    string logpath = conf.getString("logpath","./logsSwiftHohenberg");
    int numSectors = conf.getInt("numsectors",12);
    int red = conf.getInt("red",100);
    int green = conf.getInt("green",100);
    int fov = conf.getInt("fov",45);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    bool overwrite_logs = conf.getBool("overwrite_logs",1);

    //need to fix this to conform with other region files
    unsigned int numpoints;
    numpoints = 10;

    //set up a vtxVisual pointer
    vtxVisual* cv;

// initialise DRegion class setting scale
    cout << "before Dregion" << endl;
    DRegion M(scale,xspan,logpath,numpoints); //create tessellation
    cout << "after Dregion" << endl;
    M.setCreg(); //set counts to identify inner boundaries
    M.setInternalBoundary(); //set internal boundaries
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
        return -1;
    }
    else {
        cout << "Success: setInnerRegion no of outer regions " << inReg << endl;
    }
    cout << "after first setRadialSegments " << endl;
// include the analysis methods
    Analysis L;

    /*
     * now we integrated over a polygonal tesselation but this time with the equations solved
     * in each region with a separate schSolver
     */
    vector<shSolver> S;
    /*
     * Set boundaries based on polygons derived from the Dregion tesselation.
     * stored in curvedBoundary
     */
    M.populateBoundPolygon(1);
    cout << "just after setting polygonal boundaries " << M.curvedBoundary.size()<<endl;
    for (unsigned int j = 0;j<numpoints;j++) {
        S.push_back(shSolver(scale, xspan, logpath, M.curvedBoundary[j], M.centroids[j]));
        cout << "in the loop populating the schVector morph0 "<< j <<endl;
    }
    /*now set the parameters for each solver
    for (int j=0; j<numpoints;j++) {
        S[j].setParams(D_A, D_B, k1, k2, k3, k4);
    }
    */
    cout << "first setting of centroids" << endl;
    for (unsigned int j=0; j<numpoints;j++){
        std::cout  << "centroid region " << j << " is ( " << M.centroids[j].first << " , " << M.centroids[j].second << " )" << endl;
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
        std::cout  << " Region size " << M.regionHex[j].size() << endl;
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
    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
    unsigned int seed = time(NULL);
    // A rando2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);

 #ifdef COMPILE_PLOTTING
// now draw the intial tesselation
    FLT hexWidth = M.Hgrid->hexen.begin()->d/2.0;
    cerr << "d/2: " << hexWidth << endl;
    // Should the plots be saved as png images?
    const bool saveplots = conf.getBool ("saveplots", true);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    const bool vidframes = conf.getBool ("vidframes", false);
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
    cv->addLabel("Tessellation morph 0", {0.0f, 1.1f, 0.0f});
    //v1->addVisualModel(cv);
    cout << "before rendering v1 " << endl;
    v1->render();
    cout << "after rendering v1 " << endl;
    std::stringstream frame;
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    //frame << framenum++;
    frame << ".png";
    cout << " before save image " << endl;
    savePngs (logpath, "tessellation0", 0, *v1);


#endif

    // set random i.c.s
    for (unsigned int j=0; j<numpoints; j++) {
        for (auto &h : S[j].Hgrid->hexen) {
            FLT choice = ruf.get();
            S[j].psi[h.vi] = std::polar (1.0f, choice * 2.0f * 3.1415927f);
            FLT choice1 = ruf.get();
            choice1 = ruf.get();
            S[j].phi[h.vi] = std::polar (1.0f, choice1 * 2.0f * 3.1415927f);
        }
    }


 #ifdef COMPILE_PLOTTING
    // Spatial offset, for positioning of HexGridVisuals
    Vector<float> spatOff;
    float xzero = 0.0f;

    // A. Offset in x direction to the left.
    xzero = 0.7;
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    Scale<FLT,float> zscale; zscale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    Scale<FLT,float> cscale; cscale.do_autoscale = true;

    unsigned int Agrid[numpoints];

    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
            vector<FLT> regionA;
            //normalise over the region t
            regionA = L.getAbs(S[j].psi);
            Agrid[j] = v1->addVisualModel (new HexGridVisual<FLT> (v1->shaderprog,
                                                                    v1->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionA),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::Jet));
    }//end of loop over regions
    // A. Offset in x direction to the left.
    xzero -= 1.2 * M.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    unsigned int Bgrid[numpoints];
    for (unsigned int j = 0;j<numpoints;j++) { //loop over regions
            vector<FLT> regionB;
            //normalise over the region t
            regionB = L.getArgPrincipal(S[j].psi);
            Bgrid[j] = v1->addVisualModel (new HexGridVisual<FLT> (v1->shaderprog,
                                                                    v1->tshaderprog,
                                                                    S[j].Hgrid,
                                                                    spatOff,
                                                                    &(regionB),
                                                                    zscale,
                                                                    cscale,
                                                                    morph::ColourMapType::Jet));
    }//end of loop over regions
    cout << "after setting A and B grids morph0" << std::endl;
#endif
   std::cout << "after setting up graphics" << std::endl;
    // begin time stepping loop unmorphed grid solved via schSolver
    /* set up vectors for determining convergence
    std::vector<FLT> Adiff;
    std::vector<std::vector<FLT>> Apre;
    std::vector<std::vector<FLT>> Acurr;
    Adiff.resize(numpoints);
    Apre.resize(numpoints);
    Acurr.resize(numpoints);
    FLT AdiffSum = 0.0f;
  //initilise all Apre vectors to above possible field
    for (int j=0; j<numpoints;j++) {
        Apre[j].resize(S[j].A.size(),10000.0);
    }
    */
    vector<vector<complex<FLT>>> oldPsi;
    oldPsi.resize(numpoints);
    for (unsigned int j=0; j<numpoints; j++) {
        oldPsi[j].resize(S[j].n,0.0);
        oldPsi[j] = S[j].psi;
    }
    for (int i=0;i<numsteps;i++) {
        for (unsigned int j=0; j<numpoints; j++) {
            //std::cerr << "before step " << i << std::endl;
            S[j].step(dt, epsilon, g, oldPsi[j]);
            //std::cerr << "after step " << i << std::endl;
            oldPsi[j] = S[j].psi;
        } // end of time stepping
#ifdef COMPILE_PLOTTING
        if ((i % numprint) == 0) {
            for(unsigned int j=0; j<numpoints; j++) {
                vector<FLT> A;
                A.resize(S[j].n);
                A = L.getAbs(S[j].psi);
                vector<FLT> B;
                B.resize(S[j].n);
                B = L.getArgPrincipal(S[j].psi);
                FLT epsilon = 0.01;
                cerr << "max arg of normalpsi  region " << j << " is " << L.maxVal(B) << " min arg of normalpsi " << L.minVal(B) <<  " iteration " << i <<endl;
                cerr << "max val of abs(psi) region  " << j << " is " << L.maxVal(A) << " min val of abs(psi) " << L.minVal(A) <<  " iteration " << i <<endl;

                VisualDataModel<FLT>* hgv1 = (VisualDataModel<FLT>*)v1->getVisualModel (Agrid[j]);
                hgv1->updateData (&A);
                hgv1->clearAutoscaleColour();

                VisualDataModel<FLT>* hgv2 = (VisualDataModel<FLT>*)v1->getVisualModel (Bgrid[j]);
                hgv2->updateData (&B);
                hgv2->clearAutoscaleColour();
                if (saveplots) {
                    if (vidframes) {
                        savePngs (logpath, "psi", framecount, *v1);
                        ++framecount;
                    }
                    else {
                        savePngs (logpath, "psi", i , *v1);
                    }
                }

                cout << " just before isPinWheel psi complexZero" << endl;
                vector<FLT> isPinWheel = S[j].complexZero(L.getAbs(S[j].psi));
                for (unsigned int i=0; i<isPinWheel.size();i++) {
                    std::cout << "isPinwheel psi " << isPinWheel[i] << std::endl;
                }

                cout << " just after is Pinwheel psi complexZero" << endl;
                cout << " just before isPinWheel phi complexZero" << endl;
                isPinWheel.resize(0);
                isPinWheel = S[j].complexZero(L.getAbs(S[j].phi));
                for (unsigned int i=0; i<isPinWheel.size();i++) {
                    std::cout << "isPinwheel psi " << isPinWheel[i] << std::endl;
                }
                cout << " just after is Pinwheel phi complexZero" << endl;
            }//end of plotting over regions
            v1->render();
        } // end of print on plotevery
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
        steady_clock::duration sincerender = steady_clock::now() - lastrender;
        if (duration_cast<milliseconds>(sincerender).count() > 17000) { // 17 is about 60 Hz
            glfwPollEvents();
            v1->render();
            lastrender = steady_clock::now();
        }
#endif
     } //end of numsteps loop
//cout << " just after time step i = " << i << endl;

//code run at end of timestepping
    for (unsigned int j=0;j<numpoints;j++) {
        vector<FLT> B;
        B.resize(S[j].n);
        B = L.getArgPrincipal(S[j].psi);
        M.A[j] = B; //write the A fields to the DRegion array
    }
    FLT avAbsCorrelation = 0;
    cout << "just after renewcorrelate_edges morph1 " << endl;
    avAbsCorrelation = M.correlate_edges(0);
    const int max_comp = numpoints*3;
    M.random_correlate(max_comp, 0);
    cout << "just after randomcorrelate_edges morph1 " << endl;
/*first save the  ofstream outFile;
    string fname = logpath + "/first.h5";
    morph::HdfData data(fname);
    data.add_contained_vals("c",S[j].phi); //Laplacian of psi
    data.add_contained_vals("n",S[j].psi); //main variable
    data.add_val ("/g", g); //g parameter
    data.add_val ("/epsilon", epsilon); //Dn paramater
    cout << " just after writing data "  << endl;
*/
#ifdef COMPILE_PLOTTING
    cout << "Ctrl-c or press x in graphics window to exit.\n";
    v1->keepOpen();
#endif

    return 0;
};
