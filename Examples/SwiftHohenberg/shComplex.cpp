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
#include <morph/GraphVisual.h>
#include <morph/Scale.h>
#endif
#include <morph/Vector.h>
#include "shPSolver.h"
#include "analysis.h"
#include "../topo/gcal.h"
#include "../topo/analysis.h"
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
#include <chrono>
#include <complex>

#ifdef COMPILE_PLOTTING
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
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HexGrid;
using morph::HdfData;
using morph::Tools;
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
    int nonLocal = conf.getInt("nonLocal",95);
    string logpath = conf.getString("logpath","./logsSwiftHohenberg");
    int numSectors = conf.getInt("numsectors",12);
    int red = conf.getInt("red",100);
    int green = conf.getInt("green",100);
    int fov = conf.getInt("fov",45);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    bool overwrite_logs = conf.getBool("overwrite_logs",1);


// include the analysis methods
    int nbins = 50;
    int gaussBlur = -1;
    FLT radius = 1.016;
    FLT ROIwid =  0.7f*radius;
    Analysis L(nbins, gaussBlur, ROIwid);

    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
    unsigned int seed = time(NULL);
    // A rando2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);
#ifdef COMPILE_PLOTTING
    // Parameters from the config that apply only to plotting:
    const bool saveplots = conf.getBool ("saveplots", true);
    // If true, then write out the logs in consecutive order numbers,
    // rather than numbers that relate to the simulation timestep.
    const bool vidframes = conf.getBool ("vidframes", true);
    unsigned int framecount = 0;

    // Window width and height
    const unsigned int win_width = conf.getUInt ("win_width", 1025UL);
    unsigned int win_height_default = static_cast<unsigned int>(wratio * (float)win_width);
    const unsigned int win_height = conf.getUInt ("win_height", win_height_default);
    cout << "just before new Visual object" << endl;
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    morph::Visual v1(win_width, win_height, "psi and phi fields ");
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    v1.backgroundWhite();

    /*orthographic projection
    v1->ptype = morph::perspective_type::orthographic;
    v1->ortho_bl = {-1.0f, -1.0f};
    v1->ortho_tr = {1.0f, 1.0f};
    */
    // You can tweak the near and far clipping planes
    v1.zNear = 0.001;
    v1.zFar = 10000;
    // And the field of view of the visual scene.
    v1.fov = fov;
    // You can lock movement of the scene
    v1.sceneLocked = conf.getBool ("sceneLocked", false);
    // You can set the default scene x/y/z offsets
    v1.setZDefault (conf.getFloat ("z_default", -5.0f));
    v1.setSceneTransXY (x_default, y_default);
    // Make this larger to "scroll in and out of the image" faster
    v1.scenetrans_stepsize = 0.5;
    cout << "end of setting new Visual object" << endl;
    //make it current
    v1.setCurrent();

    // if using plotting, then set up the render clock
    steady_clock::time_point lastrender = steady_clock::now();

#endif
    FLT lengthScale = 29.0f;
// section for solving on the circle Tessllation
    cout << "just before creating  Solver S" << endl;
    //Readjust Dn for a single region
    pair<FLT,FLT> centroid(0.0,0.0);

/* section for solving on the circle Tessllation
// if (skipMorph) return 0;
    cout << "just before creating shCSolver" <<endl;
    shSolver S(scale, xspan, logpath, radius, centroid, lengthScale);
*/
//section for rectangular grid
    float xwidth = 2.0f;
    float ywidth = 2.0f;
    shSolver S(scale, xspan, logpath, xwidth, ywidth , lengthScale);
    int nx = std::floor(xwidth/S.ds);
    int ny = std::floor(ywidth/S.ds);
// Constructor for parallelogram domain
//    shSolver S(scale,  xspan, logpath);
    cout << "just after creating shCSolver" <<endl;
// Constructor for a rectangular domain
//    shSolver S(scale, xspan, logpath, red, green);
//    cout << "just after creating shCSolver" <<endl;
    //S.setNoFlux();
    cout << "just after setHexType()" << endl;

    string fname = logpath + "/first.h5";
//Now set i.c.s
    if (Lcontinue) {
        morph::HdfData data(fname);
        data.read_contained_vals("c",S.phi); //Laplacian of psi
        data.read_contained_vals("n",S.psi); //main variable
    }
    else {
    //random i.c.s
        for (auto &h : S.Hgrid->hexen) {
            FLT choice = ruf.get();
            S.psi[h.vi] = std::polar (1.0f, - choice * 2.0f * 3.1415927f);
            FLT choice1 = ruf.get();
            choice1 = ruf.get();
            S.phi[h.vi] = std::polar (1.0f, - choice1 * 2.0f * 3.1415927f);
        }
    }
    std::cout << "after setting ics" << std::endl;
    //
    /*set up sinusoidal initial conditions
    int di = 0;
    for (int i = 0; i<S.Hgrid->d_numrows; i++) {
        FLT eta = i * 8.0 * PI/(S.Hgrid->d_numrows);
        for (int j=0; j<S.Hgrid->d_rowlen; j++) {
            FLT theta = j * 8.0 * PI/(S.Hgrid->d_rowlen * 1.0);
            di = i*S.Hgrid->d_rowlen + j;
            S.psi[di] = std::polar (1.0f, - cos(eta) * sin(theta));
            S.phi[di] = std::polar (1.0f, - cos(eta) * sin(theta));
        }
    }
    */

    vector<FLT> psiphase;
    vector <FLT> psir;
    vector <FLT> phir;
    vector <FLT> phiphase;

    psiphase = L.getArgPrincipal(S.psi);
    psir = L.getAbs(S.psi);
    phir = L.getAbs(S.phi);
    phiphase = L.getArgPrincipal(S.phi);
    std::cout << " just before setting graphics" << std::endl;
#ifdef COMPILE_PLOTTING
    // Spatial offset, for positioning of HexGridVisuals
    morph::Vector<float> spatOff;

    // A. Offset in x direction to the left.
    // by half a hexGrid width
    float xzero = -0.5*S.Hgrid->width();
    float yzero = 0.5*S.Hgrid->width();
    float txtoff = -0.55f;

    spatOff = { xzero, yzero, 0.0 };
    morph::ColourMapType cmt = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap", "Jet"));
    // Z position scaling - how hilly/bumpy the visual will be.
    morph::Scale<FLT,float> zscale; zscale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    morph::Scale<FLT,float> cscale; cscale.do_autoscale = true;
    morph::HexGridVisual<FLT>* hgv1 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv1->setScalarData (&psiphase);
    // Z position scaling - how hilly/bumpy the visual will be.
    hgv1->zScale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    hgv1->colourScale.do_autoscale = true;
    hgv1->cm.setType (cmt);
    hgv1->hexVisMode = morph::HexVisMode::Triangles;
    hgv1->addLabel ("Psi phase", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    hgv1->finalize();
    v1.addVisualModel (hgv1);

    // A. Offset in x direction to the right.
    // move back a whole hexGrid width
    xzero += S.Hgrid->width();
    spatOff = { xzero, yzero, 0.0 };
    morph::HexGridVisual<FLT>* hgv2 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv2->setScalarData (&psir);
    // Z position scaling - how hilly/bumpy the visual will be.
    hgv2->zScale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    hgv2->colourScale.do_autoscale = true;
    hgv2->cm.setType (cmt);
    hgv2->hexVisMode = morph::HexVisMode::Triangles;
    hgv2->addLabel ("Psi modulus", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    //hgv1->addLabel ("Variable A", { -0.2f, RD.ellipse_b*-1.4f, 0.01f },
    //                morph::colour::white, morph::VisualFont::Vera, 0.1f, 48);
    // "finalize" is required before adding the HexGridVisual to the morph::Visual.
    hgv2->finalize();
    v1.addVisualModel (hgv2);

    // A. Offset in y direction down.
    // move down a whole hexGrid width
    yzero -= S.Hgrid->width();
    spatOff = { xzero, yzero, 0.0 };
    morph::HexGridVisual<FLT>* hgv3 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv3->setScalarData (&phir);
    // Z position scaling - how hilly/bumpy the visual will be.
    hgv3->zScale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    hgv3->colourScale.do_autoscale = true;
    hgv3->cm.setType (cmt);
    hgv3->hexVisMode = morph::HexVisMode::Triangles;
    hgv3->addLabel ("Phi modulus", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    hgv3->finalize();
    v1.addVisualModel (hgv3);

    // A. Offset in x direction left.
    // move back a whole hexGrid width
    xzero -= S.Hgrid->width();
    spatOff = { xzero, yzero, 0.0 };
    morph::HexGridVisual<FLT>* hgv4 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv4->setScalarData (&phiphase);
    // Z position scaling - how hilly/bumpy the visual will be.
    hgv4->zScale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    hgv4->colourScale.do_autoscale = true;
    hgv4->cm.setType (cmt);
    hgv4->hexVisMode = morph::HexVisMode::Triangles;
    hgv4->addLabel ("Phi phase", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    hgv4->finalize();
    v1.addVisualModel (hgv4);

#endif
    std::cout << "after setting up graphics" << std::endl;
    vector<complex<FLT>> oldPsi;
    oldPsi.resize(S.n, 0.0);
    S.setNonLocalR();
    std::cout << "after setting NonLocalR" << std::endl;
    S.setNonLocalC();
    std::cout << "after setting NonLocalC" << std::endl;
    oldPsi = S.psi;
    vector<bool> isPinWheel;
    int pinCount = 0;
    FLT pinDensity = 0.0;
    for (int i=0;i<numsteps;i++) {
        std::cerr << "step " << i << std::endl;
        S.step(dt, epsilon, g, oldPsi);
        std::cerr << "step " << i << std::endl;
        oldPsi = S.psi;
        FLT psiMax = 0.0;
        FLT psiMin = 0.0;
        FLT phiMax = 0.0;
        FLT phiMin = 0.0;
        if ((i % nonLocal == 0) && i>0) {
      //      S.setNonLocalR();
      //      S.setNonLocalC();
            isPinWheel.resize(0);
            bool showfft = true;
            psiphase = L.getArgPrincipal(S.psi);
            psir = L.getAbs(S.psi);
            phiphase = L.getArgPrincipal(S.phi);
            phir = L.getAbs(S.phi);
            psiMax = L.maxVal(psir);
            psiMin = L.minVal(psir);
            phiMax = L.maxVal(phir);
            phiMin = L.minVal(phir);
            cv::Mat I = L.vect2mat(psiphase, nx, ny);
            //L.getPatternFrequency(I, showfft);
            FLT low = 0.1*(psiMax-psiMin) + psiMin;
            cout << " just before isPinWheel psi complexZero" << endl;
            isPinWheel = S.complexZero(L.getAbs(S.psi),low);
            pinCount = L.countBool(isPinWheel);
            std::cout << " pinwheel count via bools " << pinCount << std::endl;
            cout << " just after is Pinwheel psi complexZero" << endl;
            cout << " just before isPinWheel phi complexZero" << endl;
            isPinWheel.resize(0);
            low = 0.1*(phiMax-psiMin) + phiMin;
            isPinWheel = S.complexZero(L.getAbs(S.phi), low);
            cout << " just after is Pinwheel phi complexZero" << endl;
            //pinDensity = pinCount/(L.columnSpacing*L.columnSpacing);
            //cout<<"pattern frequency " << L.patternFrequency << "  columnSpacing " << L.columnSpacing << " pinWheel density " << pinDensity << std::endl;
        }
#ifdef COMPILE_PLOTTING
        if ((i % numprint) == 0) {
            psiphase = L.getArgPrincipal(S.psi);
            psir = L.getAbs(S.psi);
            phiphase = L.getArgPrincipal(S.phi);
            phir = L.getAbs(S.phi);
            psiMax = L.maxVal(psir);
            psiMin = L.minVal(psir);
            phiMax = L.maxVal(phir);
            phiMin = L.minVal(phir);
            cerr << "max arg of normalpsi  " << L.maxVal(psiphase) << " min arg of normalpsi " << L.minVal(psiphase) <<  " iteration " << i <<endl;
            cout << "max arg of normalpsi  " << L.maxVal(psiphase) << " min arg of normalpsi " << L.minVal(psiphase) <<  " iteration " << i <<endl;
            cerr << "max val of abs(psi)  " << psiMax << " min val of abs(psi) " << psiMin <<  " iteration " << i <<endl;
            cout << "max val of abs(psi)  " << psiMax << " min val of abs(psi) " << psiMin <<  " iteration " << i <<endl;
            hgv1->updateData (&psiphase);
            hgv1->clearAutoscaleColour();

            hgv2->updateData (&psir);
            hgv2->clearAutoscaleColour();

            hgv3->updateData (&phir);
            hgv3->clearAutoscaleColour();

            hgv4->updateData (&phiphase);
            hgv4->clearAutoscaleColour();

            if (saveplots) {
                if (vidframes) {
                    savePngs (logpath, "psi", framecount, v1);
                    ++framecount;
                }
                else {
                    savePngs (logpath, "psi", i , v1);
                }
            }

        } // end of print on numprin     disp.closeDisplay();
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
        steady_clock::duration sincerender = steady_clock::now() - lastrender;
        if (duration_cast<milliseconds>(sincerender).count() > 170) { // 17 is about 60 Hz
            glfwPollEvents();
            v1.render();
            lastrender = steady_clock::now();
        }
#endif
     } //end of numsteps loop
//cout << " just after time step i = " << i << endl;

//code run at end of timestepping
/*
//first save the  ofstream outFile;
    morph::HdfData data(fname);
    data.add_contained_vals("c",S.phi); //Laplacian of psi
    data.add_contained_vals("n",S.psi); //main variable
    data.add_val ("/g", g); //g parameter
    data.add_val ("/epsilon", epsilon); //epsilon paramater
    cout << " just after writing data "  << endl;
//
*/
#ifdef COMPILE_PLOTTING
    cout << "Ctrl-c or press x in graphics window to exit.\n";
    v1.keepOpen();
#endif
/*
    cout << " just before isPinWheel psi" << endl;
    vector<bool> isPinWheel = S.complexZero(L.getAbs(S.psi));
    cout << " just after is Pinwheel psi " << endl;
    for (unsigned int i=0; i<isPinWheel.size();i++) {
        std::cout << "isPinwheel psi " << isPinWheel[i] << std::endl;
    }
//    for (unsigned int i=0; i<isPinWheel.size(); i++) {
//        cout << " bool value for i " << i << " is " << isPinWheel[i] << endl;
//    }

    cout << " just before isPinWheel phi" << endl;
    isPinWheel.resize(0);
    isPinWheel = S.complexZero(L.getAbs(S.phi));
    cout << " just after is Pinwheel phi " << endl;
    for (unsigned int i=0; i<isPinWheel.size();i++) {
        std::cout << "isPinwheel " << isPinWheel[i] << std::endl;
    }
    //do fft
    //std::vector<float> fitCoeffs = analysis.updateIsoORfrequencyEstimate(showFFT);

*/
    return 0;
};
