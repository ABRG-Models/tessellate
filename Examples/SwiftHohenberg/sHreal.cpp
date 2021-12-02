/*!
 * If COMPILE_PLOTTING is defined at compile time, then include the display and
 * plotting code. I usually put all the plotting code inside #defines like this so
 * that I can compile a version of the binary without plotting, for parameter searches
 * in which I am only going to be saving out HDF5 data.
 */
#ifdef COMPILE_PLOTTING
#include "region.h"
#include <morph/Visual.h>
#include <morph/HexGridVisual.h>
#include <morph/ColourMap.h>
#include <morph/VisualDataModel.h>
#include <morph/Scale.h>
#endif

#include "analysis.h"
#include "shSolver.h"
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
#include <chrono>

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
#else
        double dt = conf.getDouble("dt",0.0001);
        double epsilon = conf.getDouble("epsilon",5.0);
        double g = conf.getDouble("g",5.0);
        double xspan = conf.getDouble("xspan",5.0);
        double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
        double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
        double nnInitialOffset = conf.getDouble("nnInitialOffset",1.0);
        double ccInitialOffset = conf.getDouble("ccInitialOffset", 2.5);
#endif
    int scale = conf.getInt("scale",8);
    int numsteps = conf.getInt("numsteps",100);
    int numAdjust = conf.getInt("numAdjust",1000000);
    int numprint = conf.getInt("numprint",95);
    string logpath = conf.getString("logpath","./logsSwiftHohenberg");
    int numSectors = conf.getInt("numsectors",12);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    bool overwrite_logs = conf.getBool("overwrite_logs",1);


// include the analysis methods
    Analysis L;

    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
    unsigned int seed = time(NULL);
    // A rando2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);
#ifdef COMPILE_PLOTTING
    // Parameters from the config that apply only to plotting:
    const unsigned int plotevery = conf.getUInt ("plotevery", 10);
    // Should the plots be saved as png images?
    const bool saveplots = conf.getBool ("saveplots", true);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    const bool vidframes = conf.getBool ("vidframes", false);
    unsigned int framecount = 0;

    // Window width and height
    const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
    unsigned int win_height_default = static_cast<unsigned int>(0.8824f * (float)win_width);
    const unsigned int win_height = conf.getUInt ("win_height", win_height_default);
    cout << "just before new Visual object" << endl;
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    morph::Visual * v1;
    v1 = new morph::Visual(win_width, win_height, "psi and phi fields ");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v1->backgroundWhite();

    /*orthographic projection
    v1->ptype = morph::perspective_type::orthographic;
    v1->ortho_bl = {-1.0f, -1.0f};
    v1->ortho_tr = {1.0f, 1.0f};
    */
    // You can tweak the near and far clipping planes
    v1->zNear = 0.001;
    v1->zFar = 10000;
    // And the field of view of the visual scene.
    v1->fov = 45;
    // You can lock movement of the scene
    v1->sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v1->setZDefault (conf.getFloat ("z_default", -5.0f));
    v1->setSceneTransXY (conf.getFloat ("x_default", 0.0f),
                        conf.getFloat ("y_default", 0.0f));
    // Make this larger to "scroll in and out of the image" faster
    v1->scenetrans_stepsize = 0.5;
    cout << "end of setting new Visual object" << endl;
    //make it current
    v1->setCurrent();

    // if using plotting, then set up the render clock
    steady_clock::time_point lastrender = steady_clock::now();

#endif
// section for solving on the circle Tessllation
    cout << "just before creating circle Tessellation Solver S" << endl;
    //Readjust Dn for a single region
    FLT radius = 1.0;
    pair<FLT,FLT> centroid(0.0,0.0);

// section for solving on the circle Tessllation
// if (skipMorph) return 0;
    cout << "just before creating circle Tessellation" <<endl;
    shSolver S(scale, xspan, logpath, radius, centroid);
    FLT choice = ruf.get();
    for (auto &h : S.Hgrid->hexen) {
        if (choice > 0.5)
        {
            S.psi[h.vi] = - ruf.get() * aNoiseGain +nnInitialOffset;
            S.phi[h.vi] = - ruf.get() * aNoiseGain + ccInitialOffset;
        }
        else
        {
            S.psi[h.vi] =  ruf.get() * aNoiseGain +nnInitialOffset;
            S.phi[h.vi] =  ruf.get() * aNoiseGain + ccInitialOffset;
        }
    }
#ifdef COMPILE_PLOTTING
    // Spatial offset, for positioning of HexGridVisuals
    Vector<float> spatOff;
    float xzero = 0.0f;

    // A. Offset in x direction to the left.
    // by half a hexGrid width
    xzero -= 0.5*S.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    // Z position scaling - how hilly/bumpy the visual will be.
    Scale<FLT,float> zscale; zscale.setParams (0.2f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    Scale<FLT,float> cscale; cscale.do_autoscale = true;
    vector<FLT> normalpsi;
    normalpsi = L.normalise(S.phi);
    unsigned int Agrid = v1->addVisualModel (new HexGridVisual<FLT> (v1->shaderprog,
                                                               v1->tshaderprog,
                                                               S.Hgrid,
                                                               spatOff,
                                                               &(normalpsi),
                                                               zscale,
                                                               cscale,
                                                               morph::ColourMapType::Jet));
    //v1->getVisualModel(Agrid)->addLabel ("psi field");
    // A. Offset in x direction to the right.
    // move back a whole hexGrid width
    xzero += S.Hgrid->width();
    spatOff = { xzero, 0.0, 0.0 };
    vector<FLT> normalphi;
    normalphi = L.normalise(S.phi);
    unsigned int Bgrid = v1->addVisualModel (new HexGridVisual<FLT> (v1->shaderprog,
                                                               v1->tshaderprog,
                                                               S.Hgrid,
                                                               spatOff,
                                                               &(normalphi),
                                                               zscale,
                                                               cscale,
                                                               morph::ColourMapType::Jet));
   // v1->getVisualModel(Bgrid)->addLabel ("phi field");
#endif
    vector<FLT> oldPsi;
    oldPsi.resize(S.n, 0.0);
    oldPsi = S.psi;
    for (int i=0;i<numsteps;i++) {
        S.step(dt, epsilon, g, oldPsi);
        oldPsi = S.psi;
#ifdef COMPILE_PLOTTING
        if ((i % plotevery) == 0) {
            vector<FLT> normalpsi;
	    normalpsi.resize(S.n);
            normalpsi = L.normalise(S.psi);
            VisualDataModel<FLT>* avm = static_cast<VisualDataModel<FLT>*>(v1->getVisualModel (Agrid));
            avm->updateData (&normalpsi);
            avm->clearAutoscaleColour();

            vector<FLT> normalphi;
	    normalphi.resize(S.n);
            normalphi = L.normalise(S.phi);
            VisualDataModel<FLT>* bvm = static_cast<VisualDataModel<FLT>*>(v1->getVisualModel (Bgrid));
            bvm->updateData (&normalphi);
            bvm->clearAutoscaleColour();
            if (saveplots) {
                if (vidframes) {
                    savePngs (logpath, "nn0", framecount, *v1);
                    ++framecount;
                }
                else {
                    savePngs (logpath, "nn0", i , *v1);
                }
            }
        } // end of print on numprin     disp.closeDisplay();
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
//first save the  ofstream outFile;
    string fname = logpath + "/first.h5";
    morph::HdfData data(fname);
    data.add_contained_vals("c",S.phi); //Laplacian of psi
    data.add_contained_vals("n",S.psi); //main variable
    data.add_val ("/g", g); //g parameter
    data.add_val ("/epsilon", epsilon); //Dn paramater
    cout << " just after writing data "  << endl;

#ifdef COMPILE_PLOTTING
    cout << "Ctrl-c or press x in graphics window to exit.\n";
    v1->keepOpen();
#endif
    return 0;
};
