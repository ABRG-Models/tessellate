/*!
 * If COMPILE_PLOTTING is defined at compile time, then include the display and
 * plotting code. I usually put all the plotting code inside #defines like this so
 * that I can compile a version of the binary without plotting, for parameter searches
 * in which I am only going to be saving out HDF5 data.
 */
#ifdef COMPILE_PLOTTING
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <morph/ColourMap.h>
#include <morph/VisualDataModel.h>
#include <morph/Scale.h>
#endif

#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/Config.h>
#include <morph/Scale.h>
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
    //  open the confgig file and read in the parameters
    morph::Config conf(jsonfile);
    if (!conf.ready) {
        cerr << "Error setting up JSON config: " << conf.emsg << endl;
    }
#ifdef SINGLE
        float dt = conf.getFloat("dt",0.0001);
        float Dn = conf.getFloat("Dn",1.0);
        float Dchi = conf.getFloat("Dchi",0.0);
        float Dc = conf.getFloat("Dc",0.3);
        float xspan = conf.getFloat("xspan",5.0);
        float boundaryFalloffDist = conf.getFloat("boundaryFalloffDist",0.0078);
        float aNoiseGain = conf.getFloat("aNoiseGain",0.1);
        float nnInitialOffset = conf.getFloat("nnInitialOffet", 1.0);
        float ccInitialOffset = conf.getFloat("ccInitialOffset",2.5);
#else
        double dt = conf.getDouble("dt",0.0001);
        double Dn = conf.getDouble("Dn",1.0);
        double Dchi = conf.getDouble("Dchi",0.0);
        double Dc = conf.getDouble("Dc",0.3);
        double xspan = conf.getDouble("xspan",5.0);
        double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
        double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
        double nnInitialOffset = conf.getDouble("nnInitialOffet", 1.0);
        double ccInitialOffset = conf.getDouble("ccInitialOffset",2.5);
#endif
    int numSectors = conf.getInt("numsectors",12);
    int scale = conf.getInt("scale",8);
    int numsteps = conf.getInt("numsteps",100);
    int numprint = conf.getInt("numprint",100);
    string logpath = conf.getString("logpath", "./logsMorph") ;
    string iter = conf.getString("iter","0");
    bool LfixedSeed = conf.getBool("LfixedSeed",0);
    bool LDn = conf.getBool("LDn",0);
    //bool overwrite_logs = conf.getBool("overwrite_logs",true);
    bool skipMorph  = conf.getBool("skipMorph",false);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    unsigned int numpoints = conf.getInt("numpoints",41);
    cout << " Lcontinue " << Lcontinue << " skipMorph " << skipMorph << endl;
    ofstream afile (logpath + "/centroids.out",ios::app);
    // adjust the number of steps according to the Dn number
    //numsteps = numsteps * floor(sqrt(36.0/Dn));
    //numprint  = numprint * floor(sqrt(36.0/Dn));
    // adjust the time step for the Dn values
    dt = dt * sqrt(Dn/36.0);
    //set up a vtxVisual pointer

    unsigned int seed;
    if (LfixedSeed) {
        seed = 1;
    }
    else {
        seed = time(NULL);
    }

    // A ra2yyndo2yym uniform generator returning real/FLTing point types
    morph::RandUniform<FLT> ruf(seed);
    ofstream gfile ( logpath + "/edges.out");
    ofstream jfile ( logpath + "/results.txt");
    ofstream degfile1 (logpath + "/degree1.data");
    ofstream degfile2 (logpath + "/degree2.data");
    ofstream degfile3 (logpath + "/degree3.data");

// include the analysis methods
    Analysis L;

    /*
     * We only have a single region
     */
    ksSolver S;
/*
 * Initialise ksSolver
 */

    morph::ReadCurves r("./polygon.svg");
    BezCurvePath<FLT> bound = r.getCorticalPath();
    std::pair<FLT,FLT> centroid = std::make_pair(0.0,0.0);
    FLT radius = 1.0;
   // S = ksSolver(scale, xspan, logpath, bound, centroid);
    S = ksSolver(scale, xspan, logpath, radius, centroid);

// initialise the fields
    string fname = logpath + "/first.h5";
    cout << "just before first data read morph 0 "<< " lcontinue " << Lcontinue <<endl;
// initialise with random field
    S.Hgrid->computeDistanceToBoundary();
    if (Lcontinue)
    {
        morph::HdfData ginput(fname,1);
        cout << "just after trying to open ../logs/first.h5" << endl;
        ginput.read_contained_vals("c",S.CC);
        ginput.read_contained_vals("n",S.NN);
    }
    else
    {
	for (auto h : S.Hgrid->hexen) {
            // boundarySigmoid. Jumps sharply (100, larger is sharper) over length
            // scale 0.05 to 1. So if distance from boundary > 0.05, noise has
            // normal value. Close to boundary, noise is less.
	    FLT choice = ruf.get();
	    if (choice > 0.5)
	    {
                S.NN[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
                S.CC[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
	    }
	    else
	    {
                S.NN[h.vi] = ruf.get() * aNoiseGain +nnInitialOffset;
                S.CC[h.vi] = ruf.get() * aNoiseGain + ccInitialOffset;
	    } //end of if on +-
            if (h.distToBoundary > -0.5) { // It's possible that distToBoundary is set to -1.0
                FLT bSig = 1.0 / ( 1.0 + exp (-100.0*(h.distToBoundary- boundaryFalloffDist)) );
                //S.NN[h.vi] = (S.NN[h.vi] - nnInitialOffset) * bSig + nnInitialOffset;
                //S.CC[h.vi] = (S.CC[h.vi] - ccInitialOffset) * bSig + ccInitialOffset;
                S.NN[h.vi] = S.NN[h.vi] * bSig;
                S.CC[h.vi] = S.CC[h.vi] * bSig;
	    } //end of if on boundary distance
        }//end of loop over HexGrid
    } //end of else on Lcontinue
     cout <<  "just after field creation first morph" << endl;

 #ifdef COMPILE_PLOTTING
// now draw the intial tesselation
    // Parameters from the config that apply only to plotting:
    const unsigned int plotevery = conf.getUInt ("plotevery", 10);
    // Should the plots be saved as png images?
    const bool saveplots = conf.getBool ("saveplots", true);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    bool vidframes = conf.getBool ("vidframes", true);
    unsigned int framecount = 0;

    // Window width and height
    const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
    unsigned int win_height_default = static_cast<unsigned int>(0.8824f * (float)win_width);
    const unsigned int win_height = conf.getUInt ("win_height", win_height_default);
    cout << "just before new Visual object" << endl;
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    morph::Visual v1(win_width, win_height, "Ermentrout Keller-Segel ");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v1.backgroundWhite();
    // You can tweak the near and far clipping planes
    v1.zNear = 0.001;
    v1.zFar = 500;
    // And the field of view of the visual scene.
    v1.fov = 45;
    // You can lock movement of the scene
    v1.sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v1.setZDefault (conf.getFloat ("z_default", -10.0f));
    //v1.setSceneTransXY (conf.getFloat ("x_default", 0.69f),
    //                  conf.getFloat ("y_default", 0.68f));
    v1.setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                      conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v1.scenetrans_stepsize = 0.5;
    //add two morph::HexGridVisuals to the morph::Visual
    //a 2D scaling to apply to the visuals
    float myscale = conf.getFloat("size_scale",1.0f);
    //a position to place the hexgrid visuals
    float float_Z = 0.0f;

    cout << "end of setting new Visual object" << endl;
    //make it current
    v1.setCurrent();

    // A z position to place the hexgrid visuals
    float _Z = 0.0f;

    morph::Vector<float, 3> spatOff = {1.0, 0.25, _Z}; // spatial offset
    // Data scaling parameters
    float _m = 0.2;
    float _c = 0.0;
    morph::Scale<FLT, float> cscale;
    cscale.setParams (_m, _c);
    cout << "end of setting new Visual object" << endl;

    // Set up a 3D map of the surface RD.n[0] using a morph::HexGridVisual
//    spatOff[0] -= 0.6 * (S.Hgrid->width());
    morph::HexGridVisual<FLT>* hgv1 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv1->setSizeScale (myscale, myscale);
    hgv1->setScalarData (&S.NN);
    // You can directly set VisualDataModel::zScale and ::colourScale:
    hgv1->zScale.setParams (_m/10.0f, _c/10.0f);
    // ...or use setters to copy one in:
    hgv1->setCScale (cscale);
    hgv1->cm.setType (morph::ColourMapType::Jet);
    hgv1->hexVisMode = morph::HexVisMode::Triangles;
    //hgv1->addLabel ("n (axon density)", {-0.6f, 2.0*S.Hgrid->width(), 0},
    //                morph::colour::white, morph::VisualFont::Vera, 0.12f, 10);
    hgv1->finalize();
    v1.addVisualModel (hgv1);
    std::cout << "end of first visual block" << std::endl;

    // Set up a 3D map of the surface RD.c[0]
    spatOff[0] -= S.Hgrid->width();
    morph::HexGridVisual<FLT>* hgv2 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv2->setSizeScale (myscale, myscale);
    hgv2->setScalarData (&S.NN);
    hgv2->zScale.setParams (_m/10.0f, _c/10.0f);
    hgv2->setCScale (cscale);
    hgv2->cm.setType (morph::ColourMapType::Jet);
    //hgv2->hexVisMode = morph::HexVisMode::HexInterp;
    hgv2->hexVisMode = morph::HexVisMode::Triangles;
    //hgv2->addLabel ("c (chemoattractant)", {-0.7f, S.Hgrid->width()/2.0f, 0},
    //                morph::colour::black, morph::VisualFont::Vera, 0.04f, 64);
    hgv2->finalize();
    v1.addVisualModel (hgv2);
    // if using plotting, then set up the render clock
    std::chrono::steady_clock::time_point lastrender = steady_clock::now();
    std::cout << "end of  second visual block" << std::endl;

#endif

    std::cout << "numsteps " << numsteps << " plotevery " << plotevery <<  std::endl;
    // begin morph0 time stepping loop
    std::pair<FLT, FLT> mm; // maxmin
    for (int i=0;i<numsteps;i++) {
        S.stepEuler(dt, Dn, Dchi, Dc);
#ifdef COMPILE_PLOTTING
        if ((i % numprint) == 0) {
            //scale NN
            mm = morph::MathAlgo::maxmin (S.CC);
            std::cout << "CC range: " << std::abs(mm.second - mm.first) << std::endl;
            hgv1->colourScale.compute_autoscale (mm.second, mm.first);
            //scale lapNN
            mm = morph::MathAlgo::maxmin (S.lapCC);
            std::cout << "lapCC min: " << mm.first << " lapCC max " << mm.second << std::endl;
            hgv2->colourScale.compute_autoscale (mm.second, mm.first);
            //update data
            hgv1->updateData (&S.CC);
            hgv2->updateData (&S.lapCC);
            //hgv1->clearAutoscaleColour();
            std::cout << std::setprecision(16) << " CC " << S.sum_CC << " lapCC " << S.sum_lapCC << " time step " << i << " plotevery " << plotevery << std::endl;
            vidframes = true;
            if (i % plotevery == 0) {
                if (vidframes) {
                    savePngs (logpath, "nn0", framecount, v1);
                    ++framecount;
                }
                else {
                    savePngs (logpath, "nn0", i , v1);
                }
            }
        }
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
        steady_clock::duration sincerender = steady_clock::now() - lastrender;
        if (duration_cast<milliseconds>(sincerender).count() > 17) { // 17 is about 60 Hz
            glfwPollEvents();
            v1.render();
            lastrender = steady_clock::now();
        }
    }
  //  cerr << "Ctrl-c or press x in graphics window to exit.\n";
  //  v1->keepOpen();
#endif
    //code run at end of timestepping
    //first save the  ofstream outFile;
    cout << "just before first data write morph 0" << endl;
    morph::HdfData fdata(fname);
    fdata.add_contained_vals("c",S.CC);
    fdata.add_contained_vals("n",S.NN);
    fdata.add_val ("/Dchi", Dchi);
    fdata.add_val ("/Dn", Dn);
    fdata.add_val ("/Dc",Dc);
    cout << " just after writing data "  << endl;
    gfile << endl << "analysis on first morphing iteration " << endl;
    cout << "end of program reached successfully!" << endl;
    return 0;
};
