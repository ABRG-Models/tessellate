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
#include <morph/GraphVisual.h>
#include <morph/PolygonVisual.h>
#include <morph/Scale.h>
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

/*
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
*/
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
namespace fs = std::filesystem;

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
    std::cerr << "after reading json file " << std::endl;
 #ifdef SINGLE
        float dt = conf.getFloat("dt",0.0001);
        float epsilon = conf.getFloat("epsilon",0.1);
        float k0 = conf.getFloat("k0",10.0);
        float  g = conf.getFloat("g",5.0);
        float  xspan = conf.getFloat("xspan",5.0);
        float  boundaryFalloffDist = conf.getFloat("boundaryFalloffDist",0.0078);
        float  aNoiseGain = conf.getFloat("aNoiseGain",0.1);
        float  nnInitialOffset = conf.getFloat("nnInitialOffset",1.0);
        float  ccInitialOffset = conf.getFloat("ccInitialOffset", 2.5);
        float  x_default = conf.getFloat("x_default",0.0);
        float  y_default = conf.getFloat("y_default",0.0);
        float  wratio = conf.getFloat("wratio", 0.844f);
        float  radius = conf.getFloat("radius", 1.0);
        float  ratio = conf.getFloat("ratio", 1.0);
        float  ROIwid = conf.getFloat("ROIwid", 1.4);
        float  lengthScale = conf.getFloat("lengthScale", 29.0f);
#else
        double dt = conf.getDouble("dt",0.0001);
        double epsilon = conf.getDouble("epsilon",0.1);
        double k0 = conf.getDouble("k0",10.0);
        double g = conf.getDouble("g",5.0);
        double xspan = conf.getDouble("xspan",5.0);
        double boundaryFalloffDist = conf.getDouble("boundaryFalloffDist",0.0078);
        double aNoiseGain = conf.getDouble("aNoiseGain",0.1);
        double nnInitialOffset = conf.getDouble("nnInitialOffset",1.0);
        double ccInitialOffset = conf.getDouble("ccInitialOffset", 2.5);
        double  x_default = conf.getDouble("x_default",0.0);
        double  y_default = conf.getDouble("y_default",0.0);
        double  wratio = conf.getDouble("wratio", 0.844f);
        double radius = conf.getDouble("radius", 1.0);
        double ratio = conf.getDouble("ratio", 1.0);
        double  ROIwid = conf.getDouble("ROIwid", 1.4);
        double  lengthScale = conf.getDouble("lengthScale", 29.0f);
#endif
    int scale = conf.getInt("scale",8);
    int numsteps = conf.getInt("numsteps",100);
    int numCheck = conf.getInt("numCheck",1000000);
    int numprint = conf.getInt("numprint",95);
    int nonLocal = conf.getInt("nonLocal",95);
    string logpath = conf.getString("logpath","./logsSwiftHohenberg");
    int numSectors = conf.getInt("numsectors",12);
    int red = conf.getInt("red",100);
    int green = conf.getInt("green",100);
    int fov = conf.getInt("fov",45);
    int nbins = conf.getInt("nbins",100);
    int gaussBlur = conf.getInt("gaussBlur", -2);
    bool Lcontinue = conf.getBool("Lcontinue",false);
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    bool overwrite_logs = conf.getBool("overwrite_logs",1);
    bool showfft = conf.getBool("showfft",true);
    bool lCircle = conf.getBool("lCircle",true);
    std::cerr << "after reading the json values " << std::endl;
// include the analysis methods
    FLT halfWidth = 0.5*ROIwid;
    Analysis L(nbins, gaussBlur, ROIwid);

    ofstream pinData (logpath + "/pinCount.data",ios::app);
    cout<< "just before first data read"<< " Lcontinue " << Lcontinue <<endl;
    unsigned int seed = time(NULL);
    // A rando2yym uniform generator returning real/floating point types
    morph::RandUniform<double> ruf(seed);
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
    cout << "just before new Visual object height " << win_height << " width " << win_width << endl;
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    morph::Visual v1(win_width, win_height, "psi and phi fields ");
    cout << "just after new Visual object" << endl;
    // Set a white background . This value has the order
    // 'RGBA', though the A(alpha) makes no difference.
    //v1.backgroundWhite();

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


// section for solving on the circle Tessllation
    cout << "just before creating  Solver S" << endl;
    //Readjust Dn for a single region
    pair<FLT,FLT> centroid(0.0,0.0);



// section for solving on the circle Tessllation
// if (skipMorph) return 0;
    cout << "just before creating shCSolver" <<endl;
    shSolver S(scale, xspan, logpath, radius, centroid, lengthScale, ratio, false);
    FLT circleArea = PI*radius*radius;
    FLT hexWidth = S.Hgrid->width();
    FLT hexDepth = S.Hgrid->depth();
    FLT roIarea = ROIwid*ROIwid;
    std::cout << "hexWidth " << hexWidth << " hexDepth " << hexDepth << " circle Area " <<  circleArea << std::endl;
    cout << "just after setting boundary conditions" << endl;

    string fname = logpath + "/first.h5";
//Now set i.c.s
    vector<FLT> psiphase;
    vector<FLT> phiphase;
    vector <FLT> psir;
    vector <FLT> phir;
    morph::HdfData ginput(fname,1);
    ginput.read_contained_vals("psiR",psir); //Laplacian of psi
    ginput.read_contained_vals("psiPhase",psiphase); //main variable
    ginput.read_contained_vals("phiR",phir); //main variable
    ginput.read_contained_vals("psiPhase",phiphase); //main variable
    S.psi = L.complexify(psir, psiphase);
    S.phi = L.complexify(phir, phiphase);
    std::cout << "size of S.psi " << S.psi.size() << " size of psiR " << psir.size() << endl;
    psiphase = L.getArgPrincipal(S.psi);
    psir = L.getAbs(S.psi);
    phir = L.getAbs(S.phi);
    phiphase = L.getArgPrincipal(S.phi);
    std::cout << "size of S.psi " << S.psi.size() << " size of psiR " << psir.size() << endl;
    std::cout << "after setting i.c.s" << std::endl;
    std::cout << " just before setting graphics hexGrid width " << hexWidth << std::endl;
    // Spatial offset, for positioning of HexGridVisuals
    morph::Vector<float> spatOff;

    // A. Offset in x direction to the left.
    // by half a hexGrid width
    float xzero = -hexWidth*0.65f;
    float yzero = 0.0;
    float txtoff = -0.55f;

    spatOff = { xzero, yzero, 0.0 };
    morph::ColourMapType cmt = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap", "Rainbow"));
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
    hgv1->hexVisMode = morph::HexVisMode::HexInterp;
    //hgv1->addLabel ("Psi phase", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    hgv1->finalize();
    std::cout << "after hgv1 " << std::endl;
    v1.addVisualModel (hgv1);
    // A. Offset in x direction to the right.
    // move back a whole hexGrid width
    xzero += hexWidth;
    spatOff = { xzero, yzero, 0.0 };
    cmt = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap", "Greyscale"));
    morph::HexGridVisual<FLT>* hgv2 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv2->setScalarData (&psir);
    // Z position scaling - how hilly/bumpy the visual will be.
    hgv2->zScale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    hgv2->colourScale.do_autoscale = true;
    hgv2->cm.setType (cmt);
    hgv2->hexVisMode = morph::HexVisMode::HexInterp;
    //hgv2->addLabel ("Psi modulus", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    //hgv1->addLabel ("Variable A", { -0.2f, RD.ellipse_b*-1.4f, 0.01f },
    //                morph::colour::white, morph::VisualFont::Vera, 0.1f, 48);
    // "finalize" is required before adding the HexGridVisual to the morph::Visual.
    hgv2->finalize();
    v1.addVisualModel (hgv2);

    vector<complex<FLT>> oldPsi;
    oldPsi.resize(S.n, 0.0);
    oldPsi = S.psi;
    std::vector<bool> isPinWheel;
    S.setConvolutionIndex();
    std::cout << "after setConvolution index" << std::endl;
    S.setNonLocalR();
    std::cout << "after setting NonLocalR" << std::endl;
    S.setNonLocalC();
    std::cout << "after setting NonLocalC" << std::endl;
    int pinCount = 0;
    for (int i=0;i<numsteps;i++) {
        //std::cerr << "step " << i << std::endl;
        S.step(dt, epsilon, g, k0, oldPsi);
        oldPsi = S.psi;
        //std::cerr << "step " << i << std::endl;
    }
    FLT psiMax = 0.0;
    FLT psiMin = 0.0;
    FLT phiMax = 0.0;
    FLT phiMin = 0.0;
    isPinWheel.resize(0);
    psiphase.resize(0);
    phiphase.resize(0);
    psir.resize(0);
    phir.resize(0);
    psiphase = L.scaleVect(L.getArgPrincipal(S.psi),1.0f);
    psir = L.getAbs(S.psi);
    phiphase = L.scaleVect(L.getArgPrincipal(S.phi),1.0f);
    phir = L.getAbs(S.phi);
    psiMax = L.maxVal(psir);
    psiMin = L.minVal(psir);
    phiMax = L.maxVal(phir);
    phiMin = L.minVal(phir);
    std::pair<FLT, FLT> roIcentre;
    roIcentre.first = 0.0f; roIcentre.second = 0.0f;
    FLT low = 0.01*(psiMax-psiMin) + psiMin;
    cout << " just before isPinWheel psi complexZero" << endl;
    isPinWheel = S.complexZero(psir,low);
    std::vector<morph::Vector<FLT,3>> pWCoords = S.pinWheelCoords(isPinWheel);
    cerr << "max arg of normalpsi  " << L.maxVal(psiphase) << " min arg of normalpsi " << L.minVal(psiphase)  <<endl;
    cout << "max arg of normalpsi  " << L.maxVal(psiphase) << " min arg of normalpsi " << L.minVal(psiphase)  <<endl;
    cerr << "max val of abs(psi)  " << psiMax << " min val of abs(psi) " << psiMin << std::endl;
    cout << " just after is Pinwheel psi complexZero" << endl;
    cout << " just before isPinWheel phi complexZero" << endl;
    isPinWheel.resize(0);
    low = 0.1*(phiMax-psiMin) + phiMin;
    isPinWheel = S.complexZero(phir, low);
    cout << " just after is Pinwheel phi complexZero" << endl;
    hgv1->updateData (&psiphase);
    hgv1->clearAutoscaleColour();
    hgv2->updateData (&psir);
    hgv2->clearAutoscaleColour();
    FLT sz = S.ds;
    morph::Vector<float, 3> offset2 = { 0.0f, 0.0f, 0.0f };
    array<float,3> cl_b = morph::ColourMap<float>::jetcolour (0.78);
    //array<float,3> cl_a = morph::ColourMap<float>::jetcolour (0.0);
    array<float,3> cl_a = {0.0f, 0.0f, 0.0f};
    std::vector<morph::Vector<FLT,3>>::iterator pW;
    for (pW = pWCoords.begin(); pW != pWCoords.end(); pW++) {
        morph::Vector<float,3> vtx = *pW;
        (*pW)[0] -= hexWidth*0.65;
        vtx += morph::Vector<float, 3>({1,0,0});
        v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW, vtx, sz/4.0, 0.002f, cl_a, 6));
        (*pW)[0] += hexWidth;
        v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW, vtx, sz*2.0, 0.002f, cl_b, 6));
    }
    /*
    pW = pWCoords.begin();
    morph::Vector<float,3> vtx = *pW;
    vtx += morph::Vector<float, 3>({1,0,0});
    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW, vtx, sz/12.0f, 0.002f, cl_b, 6));
    */
    v1.render();
    int i=0;
    if (saveplots) {
        if (vidframes) {
            savePngs (logpath, "psi", framecount, v1);
             ++framecount;
        }
        else {
            savePngs (logpath, "psi", i , v1);
        }
    }
    cout << "Ctrl-c or press x in graphics window to exit.\n";
    v1.keepOpen();
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
    return 0;
};
