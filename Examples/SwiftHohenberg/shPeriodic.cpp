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
#include <morph/PolygonVisual.h>
#endif
#include <morph/Vector.h>
#include <morph/HdfData.h>
#include "shPSolver.h"
#include "analysis.h"
#include <morph/Config.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
#include <chrono>
#include <complex>

/*
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
*/
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
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
        float  ROIwid = conf.getFloat("ROIwid", 1.4);
        float  lengthScale = conf.getFloat("lengthScale", 29.0f);
        float  sigma = conf.getFloat("sigma", 0.3f);
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
        double  ROIwid = conf.getDouble("ROIwid", 1.4);
        double  lengthScale = conf.getDouble("lengthScale", 29.0f);
        double  sigma = conf.getDouble("sigmma", 0.3f);
#endif
    int scale = conf.getInt("scale",8);
    int numsteps = conf.getInt("numsteps",100);
    int numCheck = conf.getInt("numCheck",1000000);
    int numprint = conf.getInt("numprint",95);
    int nonLocal = conf.getInt("nonLocal",1000000);
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
    bool lShowPinCentres = conf.getBool("lShowPinCentres",false);
    bool lContourPeriodic = conf.getBool("lContourPeriodic",true);
    bool lPeriodic = conf.getBool("lPeriodic", true);

    std::cout << "after reading the json values g " << g << std::endl;
// include the analysis methods
    FLT halfWidth = 0.5*ROIwid;
    Analysis L(nbins, gaussBlur, ROIwid);

    ofstream pinData (logpath + "/pinCount.data",ios::app);
    ofstream maxData (logpath + "/maxVal.data",ios::app);
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
    cout << "just before new Visual object height " << win_height << " width " << win_width << endl;
    // Set up the morph::Visual object which provides the visualization scene (and
    // a GLFW window to show it in)
    morph::Visual v1(win_width, win_height, "psi and phi fields ");
    cout << "just after new Visual object" << endl;
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

/*
//section for rectangular grid
    float xwidth = 2.0f;
    float ywidth = 2.0f;
    //Constructor for rectangular domain
    shSolver S(scale, xspan, logpath, xwidth, ywidth , lengthScale, sigma);
    int nx = std::floor(xwidth/S.ds);
    int ny = S.n / nx;
    int rowlen = S.Hgrid->d_rowlen;
    int numrows = S.Hgrid->d_numrows;
    FLT hexDepth = S.Hgrid->depth();
    FLT hexWidth = S.Hgrid->width();
    FLT roIarea = ROIwid*ROIwid;
    std::cout << "nx " << nx << " ny " << ny << " rowlen " << rowlen << " numrows " << numrows <<  std::endl;
    std::cout << "width " << hexWidth << std::endl <<  " depth " << hexDepth << std::endl;
    if (rowlen*numrows != S.n) {
        std::cerr << "error rowlen*numrows " << rowlen*numrows << " S.n " << S.n << std::endl;
        //return -1;
    }
*/

// Constructor for parallelogram domain
    std::cout << "first the values for the corners reported from morphologica methods" << std::endl;
    shSolver S(scale,  xspan, logpath, lengthScale, sigma, lPeriodic);
    FLT cos60 = morph::SQRT_OF_3_OVER_2_F;
    std::cout << "now the values from the actual parallelogram HexGrid constructed by morphologica cos60 " << cos60 << std::endl;
    FLT pspan = xspan/3.0f;
    int numrows = S.Hgrid->d_numrows;
    int rowlen =  S.Hgrid->d_rowlen;
    FLT hexWidth = S.Hgrid->width();
    FLT hexDepth = S.Hgrid->depth();
    FLT parDepth = (numrows-1) * S.ds;
    FLT cutWidth = (numrows-1) * S.ds / 2.0;
    FLT parWidth = (rowlen-1) * S.ds;
    FLT parArea = parDepth*parWidth;
    FLT hexArea = hexWidth * hexDepth;
    FLT roIarea = ROIwid*ROIwid;
    int topLeftIndex = (numrows-1)*rowlen;
    int topRightIndex = numrows*rowlen-1;
    std::cout << "parallelogram solver rowlen " << rowlen << " numrows  " << numrows << " hex spacing " << S.ds << " num of Hexes " << S.Hgrid->num() << std::endl;
    std::cout << "Hex width " << hexWidth << " Hex Depth " << hexDepth << std::endl;
    std::cout << "Par width " << parWidth << " cutWidth " << cutWidth << " Par depth " << parDepth << " hexArea " << hexArea << " parArea " << parArea << std::endl;
    std::cout << "bot left " << S.Hgrid->d_x[0] << " , " << S.Hgrid->d_y[0] << std::endl;
    std::cout << "bot right " << S.Hgrid->d_x[rowlen-1] << " , " << S.Hgrid->d_y[rowlen-1] << std::endl;
    std::cout << "top left " << S.Hgrid->d_x[topLeftIndex] << " , " << S.Hgrid->d_y[topLeftIndex] << std::endl;
    std::cout << "top right " << S.Hgrid->d_x[topRightIndex] << " , " << S.Hgrid->d_y[topRightIndex] << std::endl;
/*
    if ((numrows-1)%4 == 0) {
        S.setPeriodicEven();
    }
    else {
        S.setPeriodicOdd();
    }
*/
    cout << "just after setting boundary conditions" << endl;

    string fname = logpath + "/first.h5";
//Now set i.c.s
    vector<FLT> psiphase;
    vector<FLT> phiphase;
    vector <FLT> psir;
    vector <FLT> phir;
    cout << "just before first data read morph 0 " << endl;
    if (Lcontinue) {
        morph::HdfData ginput(fname, true);
        cout << "just before HdfData call" << endl;
        ginput.read_contained_vals("psiR",psir); //Laplacian of psi
        ginput.read_contained_vals("psiPhase",psiphase); //main variable
        ginput.read_contained_vals("phiR",phir); //main variable
        ginput.read_contained_vals("psiPhase",phiphase); //main variable
        S.psi = L.complexify(psir, psiphase);
        S.phi = L.complexify(phir, phiphase);
        std::cout << "size of S.psi " << S.psi.size() << " size of psiR " << psir.size() << endl;
    }
    else {
    //

    //random i.c.s
        for (auto &h : S.Hgrid->hexen) {
            FLT choice = ruf.get();
            S.psi[h.vi] = std::polar (1.0f, - choice * 2.0f * 3.1415927f);
            FLT choice1 = ruf.get();
            choice1 = ruf.get();
            S.phi[h.vi] = std::polar (1.0f, - choice1 * 2.0f * 3.1415927f);
        }
    }
    cout << "just after first data read morph 0 " << endl;
    //
    /*set up sinusoidal initial conditions parallelogram

    int di = 0;
    FLT etaInc = 32.0*PI/(1.0f*(numrows-1));
    FLT thetaInc = 32.0*PI/(1.0f*(rowlen-1));
        for (int i = 0; i<numrows; i++) {
            FLT eta = i * etaInc;
            for (int j=0; j<rowlen; j++) {
                FLT theta = j*thetaInc;
                S.psi[di] = std::polar (1.0f, - theta + eta);
                S.phi[di] = std::polar (1.0f, -  theta + eta);
                di++;
                //std::cout << " in i.c loop i " << i << " j " << j << " di " << di << std::endl;
            }
        }

    */
    /*
    //set up sinusoidal initial conditions rectangle
        int di = 0;
        int tcount = 0;
        for (int i = 0; i<nx; i++) {
            FLT eta = i * 2.0 * PI/(nx * 1.0);
            std::cout << "eta = " << eta << std::endl;
            for (int j=0; j<ny; j++) {
                FLT theta = j * 2.0 * PI/(ny * 1.0);
                //std::cout << "eta = " << eta << std::endl;
                di = i * nx + j;
                tcount++;
                S.psi[di] = std::polar (1.0f, - cos(eta) * sin(theta));
                S.phi[di] = std::polar (1.0f, - cos(eta) * sin(theta));
            }
        }
        std::cout << "number of hexes given psi values " << tcount << std::endl;
        }
    */

        psiphase = L.getArgPrincipal(S.psi);
        psir = L.getAbs(S.psi);
        phir = L.getAbs(S.phi);
        phiphase = L.getArgPrincipal(S.phi);
        vector<FLT> psiReal = L.getReal(S.psi);
        vector<FLT> psiImag = L.getImag(S.psi);
        std::cout << "size of S.psi " << S.psi.size() << " size of psiR " << psir.size() << endl;
    std::cout << "after setting i.c.s" << std::endl;
    std::cout << " just before setting graphics hexGrid width " << hexWidth << std::endl;
#ifdef COMPILE_PLOTTING
    // Spatial offset, for positioning of HexGridVisuals
    morph::Vector<float> spatOff;

    // A. Offset in x direction to the left.
    // by half a hexGrid width
    float xzero = -0.4*hexWidth;
    float yzero = 0.0;
    float txtoff = -0.5f;

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
//    hgv1->addLabel ("Psi phase", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    std::cout << "after hgv1 " << std::endl;
    hgv1->finalize();
    v1.addVisualModel (hgv1);
    // A. Offset in x direction to the right.
    // move back a whole hexGrid width
    xzero += hexWidth*0.8;
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
    hgv2->addLabel ("Psi modulus", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    //hgv1->addLabel ("Variable A", { -0.2f, RD.ellipse_b*-1.4f, 0.01f },
    //                morph::colour::white, morph::VisualFont::Vera, 0.1f, 48);
    // "finalize" is required before adding the HexGridVisual to the morph::Visual.
    hgv2->finalize();
    v1.addVisualModel (hgv2);

    std::cout << "after hgv2 " << std::endl;
    /*
    // A. Offset in y direction down.
    // move down a whole hexGrid width
    yzero -= hexWidth + 0.1f;
    spatOff = { xzero, yzero, 0.0 };
    cmt = morph::ColourMap<FLT>::strToColourMapType (conf.getString ("colourmap", "Rainbow"));
    morph::HexGridVisual<FLT>* hgv3 = new morph::HexGridVisual<FLT> (v1.shaderprog, v1.tshaderprog, S.Hgrid, spatOff);
    hgv3->setScalarData (&phiphase);
    // Z position scaling - how hilly/bumpy the visual will be.
    hgv3->zScale.setParams (0.0f, 0.0f);
    // The second is the colour scaling. Set this to autoscale.
    hgv3->colourScale.do_autoscale = true;
    hgv3->cm.setType (cmt);
    hgv3->hexVisMode = morph::HexVisMode::HexInterp;
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


    // Graph of frequency estimate
    std::vector<float> graphX(1,0);
    std::vector<float> graphY(1,0);
    std::vector<float> graphX2(1,0);
    std::vector<float> graphY2(1,0);
    std::vector<float> graphX3(2,0);
    std::vector<float> graphY3(2,0);
    graphY3[1] = 1.0;
    int sampwid = nbins;
    float wid = 2.0;
    float hei = 2.0;
    FLT grid2offx = xzero-1.5*hexWidth;
    morph::GraphVisual<float>* gvPinDensity = new morph::GraphVisual<float> (v1.shaderprog, v1.tshaderprog, morph::Vector<float>{grid2offx,-hexWidth*0.7,0.0f});
    morph::DatasetStyle ds;
    ds.linewidth = 0.00;
    ds.linecolour = {0.0, 0.0, 0.0};
    ds.markerstyle = morph::markerstyle::circle;
    ds.markersize = 0.02;
    ds.markercolour = {0.0, 0.0, 0.0};
    ds.markergap = 0.0;
    gvPinDensity->xlabel="frequency (cycles/ROI-width)";
    gvPinDensity->ylabel="FFT magnitude";
    gvPinDensity->setsize(wid,hei);
    gvPinDensity->setlimits (0,(float)sampwid*1.0,0,1.0); // plot up to nyquist (pixels / 2)
    gvPinDensity->setdata (graphX, graphY, ds);
    morph::DatasetStyle ds2;
    ds2.markerstyle = morph::markerstyle::circle;
    ds2.markersize = 0.0;
    ds2.markercolour = {0.0, 0.0, 0.0};
    ds2.markergap = 0.0;
    ds2.linewidth = 0.01;
    ds2.linecolour = {1.0, 0.0, 0.0};
    gvPinDensity->setdata (graphX2, graphY2, ds2);
    morph::DatasetStyle ds3;
    ds3.markerstyle = morph::markerstyle::circle;
    ds3.markersize = 0.0;
    ds3.markercolour = {0.0, 0.0, 0.0};
    ds3.markergap = 0.0;
    ds3.linewidth = 0.01;
    ds3.linecolour = {0.0, 0.0, 1.0};
    gvPinDensity->setdata (graphX3, graphY3, ds3);
    gvPinDensity->finalize();
    v1.addVisualModel (static_cast<morph::VisualModel*>(gvPinDensity));
    */
#endif
    std::cout << "after setting up graphics" << std::endl;
    vector<complex<FLT>> oldPsi;
    oldPsi.resize(S.n, 0.0);
    S.setConvolutionIndex();
    std::cout << "after setConvolution index" << std::endl;
    S.setNonLocalR();
    std::cout << "after setting NonLocalR" << std::endl;
    S.setNonLocalC();
    std::cout << "after setting NonLocalC" << std::endl;
    /*
    if ((numrows-1)%4 == 0) {
        S.setPeriodicEven();
    }
    else {
        S.setPeriodicOdd();
    }
    */
    oldPsi = S.psi;
    vector<bool> isPinWheel;
    int pinCount = 0;
    FLT pinDensity = 0.0;
    cv::Mat I;
    vector<vector<FLT>> psiImg;
    //morph::HdfData outdata(fname,morph::FileAccess::ReadWrite);
    morph::HdfData outdata(fname,false);
    unsigned int chkCount = 0;
    for (int i=0;i<numsteps;i++) {
        //std::cerr << "step " << i << std::endl;
        S.step(dt, epsilon, g, k0, oldPsi);
        //std::cerr << "step " << i << std::endl;
        oldPsi = S.psi;
        FLT psiMax = 0.0;
        FLT psiMin = 0.0;
        FLT phiMax = 0.0;
        FLT phiMin = 0.0;
        if ((i % nonLocal == 0) && i>0) {
            S.setNonLocalR();
            S.setNonLocalC();
        }
        if (i % numprint == 0) {
            isPinWheel.resize(0);
            psiphase.resize(0);
            phiphase.resize(0);
            psir.resize(0);
            phir.resize(0);
            FLT fscale = 1.0/(2.0*PI);
            psiphase = L.scaleVect(L.getArgPrincipal(S.psi),fscale);
            psir = L.getAbs(S.psi);
            phiphase = L.scaleVect(L.getArgPrincipal(S.phi),scale);
            phir = L.getAbs(S.phi);
            psiReal = L.getReal(S.psi);
            psiImag = L.getImag(S.psi);
            psiMax = L.maxVal(psir);
            psiMin = L.minVal(psir);
            phiMax = L.maxVal(phir);
            phiMin = L.minVal(phir);
            cout << "psir top left " << psir[(numrows-1)*rowlen] << " psiphase top left  " << psiphase[(numrows-1)*rowlen] << std::endl;
            cout << "phir top left  " << phir[(numrows-1)*rowlen] << " phiphase top left  " << phiphase[(numrows-1)*rowlen] << std::endl;
            cout << "psir bottom left " << psir[0] << " psiphase bottom left  " << psiphase[0] << std::endl;
            cout << "phir bottom left  " << phir[0] << " phiphase bottom left  " << phiphase[0] << std::endl;
            cout << "psir top right " << psir[numrows*rowlen-1] << " psiphase top right  " << psiphase[numrows*rowlen-1] << std::endl;
            cout << "phir top right  " << phir[numrows*rowlen-1] << " phiphase top right  " << phiphase[numrows*rowlen-1] << std::endl;
            cout << "psir bottom right " << psir[rowlen-1] << " psiphase bottom right  " << psiphase[rowlen-1] << std::endl;
            cout << "phir bottom right  " << phir[rowlen-1] << " phiphase bottom right  " << phiphase[rowlen-1] << std::endl;
            cout << "psir middle " << psir[numrows/2 + rowlen/2] << " psiphase middle  " << psiphase[numrows/2 + rowlen/2] << std::endl;
            cout << "phir middle  " << phir[numrows/2 + rowlen/2] << " phiphase middle  " << phiphase[numrows/2 + rowlen/2] << std::endl;
            std::pair<FLT, FLT> roIcentre;
            roIcentre.first = 0.0f; roIcentre.second = 0.0f;
            int iRoI = floor(ROIwid/S.ds) + 1; //rowlength of square
            //get the field limited to the square ROI
            vector<vector<FLT>> roIfield= S.fieldInRoI(psiphase, iRoI, halfWidth, roIcentre);
            //I = L.vect2matCut(psiphase, numrows, rowlen);
            //convert teh field to a mat format
            I = L.sqmatrix2mat(roIfield, iRoI);
            //here is the ROI field in vector form
            psiImg = roIfield;
            std::cout << "just after cv:Mat" << std::endl;
            //get the pattern frequency, if showfft show the full 2D spectrum
            FLT freqscale = (scale - 7)*1.0f + 1.0f;
            std::cout << "scaling " << freqscale << std::endl;
            std::vector<FLT> fitCoeffs = L.getPatternFrequency(I, showfft, freqscale);
            if (showfft) {
                cv::imshow("my window",I);
                cv::waitKey();
            }
            std::cout << "just after get pattern frequency" << std::endl;
            //for finding the isolated minima of the r component of psi
            FLT low = 0.1*(psiMax-psiMin) + psiMin;
            cout << " just before isPinWheel via intersect method" << endl;
            //get the contours of real and imag = 0
            //Periodic is laid out in row format only applied to parallelogram or rectangle
            //Spiral, hexes are build out from the centre
            vector<bool> realZero;
            vector<bool> imagZero;
            if (lContourPeriodic) {
                realZero = S.isContourZeroPeriodic(psiReal);
                imagZero = S.isContourZeroPeriodic(psiImag);
            }
            else {
                realZero = S.isContourZeroSpiral(psiReal);
                imagZero = S.isContourZeroSpiral(psiImag);
            }
            //return a bool mask array where they intersect, count them and get their coords
            //isPinWheel = S.intersectPeriodic(realZero, imagZero);
            isPinWheel = S.intersectNoFlux(realZero, imagZero);
            std::vector<morph::Vector<FLT,3>> pWCoords3 = S.pinWheelCoords(isPinWheel);
            std::vector<morph::Vector<FLT,3>> pWCoords = S.pinWheelCoords(imagZero);
            std::vector<morph::Vector<FLT,3>> pWCoords1 = S.pinWheelCoords(realZero);
            //do we count in the whole region or in an RoI
            pinCount = L.countBool(isPinWheel);
            std::cout << " pinwheel count psi intersect in whole region " << pinCount << std::endl;
            pinCount = S.pinCountInRoI(isPinWheel, halfWidth, roIcentre);
            std::cout << " pinwheel count psi intersect in RoI " << pinCount << std::endl;
            //now get the pinwheel count by minima of the r component
            isPinWheel.resize(0);
            isPinWheel = S.complexZero(psir,low);
            std::cout << "pincount via minima of psi.r " << L.countBool(isPinWheel) << std::endl;
            //these will be the coords of the minima of r
            //scaling factor L.columnSpacing is multiplied by length of RoI
            //scaling factor pinWheel count is divided by area of counting
            FLT columnSpacing = L.columnSpacing;
            pinDensity = pinCount*columnSpacing*columnSpacing/roIarea;
            pinData << " " << pinCount << " " << pinDensity << std::endl;
            cout<<"pattern frequency " << L.patternFrequency << "  columnSpacing " << columnSpacing << " pinWheel density " << pinDensity << std::endl;
            cerr<<"pattern frequency " << L.patternFrequency << "  columnSpacing " << columnSpacing << " pinWheel density " << pinDensity << std::endl;
            cerr << "max arg of normalpsi  " << L.maxVal(psiphase) << " min arg of normalpsi " << L.minVal(psiphase) <<  " iteration " << i <<endl;
            cout << "max arg of normalpsi  " << L.maxVal(psiphase) << " min arg of normalpsi " << L.minVal(psiphase) <<  " iteration " << i <<endl;
            cerr << "max val of abs(psi)  " << psiMax << " min val of abs(psi) " << psiMin <<  " iteration " << std::endl;
            maxData << psiMax << " " << psiMin  << std::endl;
            cout << " just after is Pinwheel psi complexZero" << endl;
            cout << " just before isPinWheel phi complexZero" << endl;
            isPinWheel.resize(0);
            low = 0.1*(phiMax-psiMin) + phiMin;
            isPinWheel = S.complexZero(phir, low);
            cout << " just after is Pinwheel phi complexZero" << endl;
#ifdef COMPILE_PLOTTING
            array<float,3> cl_a = {0.0f, 0.0f, 1.0f};
            array<float,3> cl_b = {1.0f, 0.0f, 0.0f};
            array<float,3> cl_c = {0.0f, 0.0f, 0.0f};
            std::vector<morph::Vector<FLT,3>>::iterator pW;
            std::vector<morph::Vector<FLT,3>>::iterator pW1;
            std::vector<morph::Vector<FLT,3>>::iterator pW3;
            morph::Vector<float, 3> offset2 = { 0.0f, 0.0f, 0.0f };
            FLT sz = S.ds/2.0f;
            int pwCount = 0;
            std::cout << "just before lShowPinCentres" << std::endl;
            if (lShowPinCentres) {
                for (pW = pWCoords.begin(); pW != pWCoords.end(); pW++) {
                    morph::Vector<float,3> vtx = *pW;
                    (*pW)[0] -= hexWidth*0.4;
                    (*pW)[1] += 0.0;
                    vtx += morph::Vector<float, 3>({1,0,0});
                    //v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW, vtx, sz, 0.01f, cl_c, 6));
                    (*pW)[0] += hexWidth*0.8;
                    (*pW)[1] += hexWidth*0.0;
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW, vtx, sz*1.5, 0.01f, cl_b, 6));
                }
                std::cout << "just after pW in  lShowPinCentres" << std::endl;
                for (pW1 = pWCoords1.begin(); pW1 != pWCoords1.end(); pW1++) {
                    morph::Vector<float,3> vtx = *pW1;
                    (*pW1)[0] -= hexWidth*0.4;
                    (*pW1)[1] += 0.0;
                    vtx += morph::Vector<float, 3>({1,0,0});
                    //v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW1, vtx, sz, 0.01f, cl_c, 6));
                    (*pW1)[0] += hexWidth*0.8;
                    (*pW1)[1] += hexWidth*0.0;
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW1, vtx, sz*1.5, 0.01f, cl_a, 6));
                }
                std::cout << "just after pWi1 in  lShowPinCentres" << std::endl;
                /*
                for (pW3 = pWCoords3.begin(); pW3 != pWCoords3.end(); pW3++) {
                    morph::Vector<float,3> vtx = *pW3;
                    (*pW3)[0] -= hexWidth*0.6;
                    (*pW3)[1] += 0.0;
                    vtx += morph::Vector<float, 3>({1,0,0});
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW3, vtx, sz, 0.01f, cl_c, 6));
                    (*pW3)[0] += hexWidth+0.10f;
                    (*pW3)[1] += hexWidth*0.0;
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW3, vtx, sz, 0.0f, cl_c, 6));
                    pwCount ++;
                }
                */

                std::cout << "pinWheel count via pWCoords3 " << pwCount << std::endl;
            }
            hgv1->updateData (&psiphase);
            hgv1->clearAutoscaleColour();

            hgv2->updateData (&psir);
            hgv2->clearAutoscaleColour();
/*
            hgv3->updateData (&phiphase);
            hgv3->clearAutoscaleColour();


            graphX = L.xs;
            graphY = L.ys;
            for (int i=0;i<nbins;i++){
                std::cout << " frequency " << graphX[i] << " power " << graphY[i] << std::endl;
            }
            int nsamp = nbins-10;
            float xmax = nsamp;

            graphX2.resize(nsamp,0.0);
            graphY2.resize(nsamp,0.0);
            arma::vec xfit(nsamp);
            for(int i=0;i<nsamp;i++){
                graphX2[i] = xmax*(float)i/(float)(nsamp-1);
                xfit[i] = graphX2[i];
            }
            arma::vec cf(fitCoeffs.size());
            for(int i=0;i<fitCoeffs.size();i++){
                cf[i] = fitCoeffs[i];
            }
            arma::vec yfit = arma::polyval(cf,xfit);
            graphY2.resize(nsamp,0);
            for(int i=0;i<nsamp;i++){
                graphY2[i] = yfit[i];
            }
            for (int i=0;i<nsamp;i++){
                std::cout << " frequency via fit " << graphX2[i] << " power " << graphY2[i] << std::endl;
            }

            graphX3[0] = L.patternFrequency;
            graphX3[1] = L.patternFrequency;

            gvPinDensity->update (graphX, graphY, 0);
            gvPinDensity->update (graphX2, graphY2, 1);
            gvPinDensity->update (graphX3, graphY3, 2);
*/
            if (saveplots) {
                if (vidframes) {
                    savePngs (logpath, "psi", framecount, v1);
                    ++framecount;
                }
                else {
                    savePngs (logpath, "psi", i , v1);
                }
            }
        // rendering the graphics. After each simulation step, check if enough time
        // has elapsed for it to be necessary to call v1.render().
            steady_clock::duration sincerender = steady_clock::now() - lastrender;
            if (duration_cast<milliseconds>(sincerender).count() > 170) { // 17 is about 60 Hz
                glfwPollEvents();
                v1.render();
                lastrender = steady_clock::now();
            }
        //
#endif
        } //end of if on numprint

//code for checkpoint
        if (i%numCheck == numCheck-1) {
            std::cout << "in checkpointing " << std::endl;
            psir.resize(0);
            psiphase.resize(0);
            phir.resize(0);
            phiphase.resize(0);
            psiphase = L.scaleVect(L.getArgPrincipal(S.psi), 1.0f);
            psir = L.getAbs(S.psi);
            phiphase = L.scaleVect(L.getArgPrincipal(S.phi),1.0f);
            phir = L.getAbs(S.phi);
            std::cout << " just before first write " << std::endl;
            string strpsiR = "psiR" + to_string(chkCount);
            string strpsiPhase = "psiPhase" + to_string(chkCount);
            string strphiR = "psiR" + to_string(chkCount);
            string strphiPhase = "psiPhase" + to_string(chkCount);
            char chrpsiR[strpsiR.size() + 1];
            strcpy(chrpsiR, strpsiR.c_str());
            char chrpsiPhase[strpsiPhase.size() + 1];
            strcpy(chrpsiPhase, strpsiPhase.c_str());
            char chrphiR[strphiR.size() + 1];
            strcpy(chrphiR, strphiR.c_str());
            char chrphiPhase[strphiPhase.size() + 1];
            strcpy(chrphiPhase, strphiPhase.c_str());
            outdata.add_contained_vals(chrpsiPhase, psiphase); //main variable
            outdata.add_contained_vals(chrphiPhase, phiphase); //Laplacian
            outdata.add_contained_vals(chrpsiR, psir); //Main variable
            outdata.add_contained_vals(chrphiR, phir); //Laplacian
    //outdata.add_contained_vals("psiImg", psiImg); //psi as a matrix
    //data.add_val ("/g", g); //g parameter
    //data.add_val ("/epsilon", epsilon); //epsilon paramater
            cout << " just after writing data chkCount "  << chkCount <<  endl;
            chkCount++;
        }

     } //end of numsteps loop
     cout << " end of timesteps loop " << endl;

//code run at end of timestepping
//first save the  ofstream outFile;
    psir.resize(0);
    psiphase.resize(0);
    phir.resize(0);
    phiphase.resize(0);
    psiphase = L.scaleVect(L.getArgPrincipal(S.psi), 1.0f);
    psir = L.getAbs(S.psi);
    phiphase = L.scaleVect(L.getArgPrincipal(S.phi),1.0f);
    phir = L.getAbs(S.phi);
    std::cout << " just before first write " << std::endl;
    outdata.add_contained_vals("psiPhase", psiphase); //main variable
    outdata.add_contained_vals("phiPhase", phiphase); //Laplacian
    outdata.add_contained_vals("psiR", psir); //Main variable
    outdata.add_contained_vals("phiR", phir); //Laplacian
    //outdata.add_contained_vals("psiImg", psiImg); //psi as a matrix
    //data.add_val ("/g", g); //g parameter
    //data.add_val ("/epsilon", epsilon); //epsilon paramater
    cout << " just after writing data end of prog chkCount "  << chkCount <<  endl;
//
/*
    cout << "Ctrl-c or press x in graphics window to exit.\n";
    v1.keepOpen();
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
