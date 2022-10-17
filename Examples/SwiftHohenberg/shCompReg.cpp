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
        float  ratio = conf.getFloat("ratio", 2.0);
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
        double radius = conf.getDouble("radius", 1.0);
        double ratio = conf.getDouble("ratio", 1.0);
        double  ROIwid = conf.getDouble("ROIwid", 1.4);
        double  lengthScale = conf.getDouble("lengthScale", 29.0f);
        double  sigma = conf.getDouble("sigma", 0.3f);
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
    bool lCircle = conf.getBool("lCircle",false);
    bool LfixedSeed = conf.getBool("LfixedSeed",false);
    bool overwrite_logs = conf.getBool("overwrite_logs",1);
    bool showfft = conf.getBool("showfft",true);
    bool lContourPeriodic = conf.getBool("lContourPeriodic",true);
    bool lShowPinCentres = conf.getBool("lShowPinCentres",false);
    bool lPeriodic = conf.getBool("lPeriodic", true);

    std::cerr << "after reading the json values " << std::endl;
// include the analysis methods
    FLT halfWidth = 0.5*ROIwid;
    Analysis L(nbins, gaussBlur, ROIwid);
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
    cout<< "just before first data read"<< " Lcont


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
// section for solving on the circle Tessllation
    cout << "just before creating  Solver S" << endl;
    //Readjust Dn for a single region
    pair<FLT,FLT> centroid(0.0,0.0);

// section for solving on the circle Tessllation
// if (skipMorph) return 0;
    cout << "just before creating shCSolver" <<endl;
    shSolver S(scale, xspan, logpath, radius, centroid, lengthScale, sigma);
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
    if (Lcontinue) {
        morph::HdfData ginput(fname,1);
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
    //

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
        psiphase = L.getArgPrincipal(S.psi);
        psir = L.getAbs(S.psi);
        phir = L.getAbs(S.phi);
        phiphase = L.getArgPrincipal(S.phi);
        std::cout << "size of S.psi " << S.psi.size() << " size of psiR " << psir.size() << endl;
    }
    std::cout << "after setting i.c.s" << std::endl;
    std::cout << " just before setting graphics hexGrid width " << hexWidth << std::endl;
#ifdef COMPILE_PLOTTING
    // Spatial offset, for positioning of HexGridVisuals
    morph::Vector<float> spatOff;

    // A. Offset in x direction to the left.
    // by half a hexGrid width
    float xzero = -0.6*hexWidth;
    float yzero = 0.6*hexWidth;
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

    // move back a whole hexGrid width
    xzero += hexWidth + 0.1f;
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
    // A. Offset in y direction down.
    // move down a whole hexGrid width
    yzero -= hexDepth + 0.1f;
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
    //hgv3->addLabel ("Phi modulus", {-0.05f, txtoff, 0.0f}, morph::colour::black, morph::VisualFont::VeraSerif, 0.05, 56);
    hgv3->finalize();
    v1.addVisualModel (hgv3);

    // A. Offset in x direction left.
    // move back a whole hexGrid width
    /*
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
    */


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
    FLT grid2offx = xzero- 1.6*hexDepth;
    morph::GraphVisual<float>* gvPinDensity = new morph::GraphVisual<float> (v1.shaderprog, v1.tshaderprog, morph::Vector<float>{grid2offx, -hexDepth, 0.0f});
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
    gvPinDensity->setlimits (0,(float)sampwid,0,1.0); // plot up to nyquist (pixels / 2)
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
    oldPsi = S.psi;
    vector<bool> isPinWheel;
    int pinCount = 0;
    FLT pinDensity = 0.0;
    cv::Mat I;
    vector<vector<FLT>> psiImg;
    //morph::HdfData outdata(fname,morph::FileAccess::ReadWrite);
    morph::HdfData outdata(fname,false);
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
            psiphase = L.scaleVect(L.getArgPrincipal(S.psi),1.0f);
            psir = L.getAbs(S.psi);
            phiphase = L.scaleVect(L.getArgPrincipal(S.phi),1.0f);
            phir = L.getAbs(S.phi);
            vector<FLT> psiReal = L.getReal(S.psi);
            vector<FLT> psiImag = L.getImag(S.psi);
            psiMax = L.maxVal(psir);
            psiMin = L.minVal(psir);
            phiMax = L.maxVal(phir);
            phiMin = L.minVal(phir);
            std::pair<FLT, FLT> roIcentre;
            roIcentre.first = 0.0f; roIcentre.second = 0.0f;
            int iRoI = 2*floor(halfWidth/S.ds) + 1; //rowlength of square
            vector<vector<FLT>> roIfield= S.fieldInRoI(psiphase, iRoI, halfWidth, roIcentre);
            //I = L.vect2matCut(psiphase, numrows, rowlen);
            I = L.sqmatrix2mat(roIfield, iRoI);
            psiImg = roIfield;
            std::cout << "just after cv:Mat" << std::endl;
            std::vector<FLT> fitCoeffs = L.getPatternFrequency(I, showfft);
            if (showfft) {
                cv::imshow("my window",I);
                cv::waitKey();
            }
            std::cout << "just after get pattern frequency " << L.patternFrequency << std::endl;
            FLT low = 0.1*(psiMax-psiMin) + psiMin;
            cout << " just before isPinWheel psi complexZero" << endl;
            //detect pinWheel by crossing of real and imaginary contours
            vector<bool> realZero;
            vector<bool> imagZero;
            if (lContourPeriodic){
                realZero = S.isContourZeroPeriodic(psiReal);
                imagZero = S.isContourZeroPeriodic(psiImag);
            }
            else {
                realZero = S.isContourZeroSpiral(psiReal);
                imagZero = S.isContourZeroSpiral(psiImag);
            }
            cout << " just after isPinWheel psi complexZero" << endl;
            isPinWheel = S.intersectNoFlux(realZero, imagZero);
            std::vector<morph::Vector<FLT,3>> pWCoords3 = S.pinWheelCoords(isPinWheel);
            std::vector<morph::Vector<FLT,3>> pWCoords = S.pinWheelCoords(imagZero);
            std::vector<morph::Vector<FLT,3>> pWCoords1 = S.pinWheelCoords(realZero);
            std::cout << "pincount psi intersects " << L.countBool(isPinWheel) << std::endl;
            pinCount = S.pinCountInRoI(isPinWheel, halfWidth, roIcentre);
            //do we count in the whole region or in an RoI
            pinCount = L.countBool(isPinWheel);
            std::cout << " pinwheel count psi intersect in whole region " << pinCount << std::endl;
            std::cout << " pinwheel count psi intersect in RoI" << pinCount << std::endl;
            //now get the pinwheel count by minima of the r component
            isPinWheel.resize(0);
            isPinWheel = S.complexZero(psir,low);
            std::cout << "pincount vim minima of psi.r " << L.countBool(isPinWheel) << std::endl;
            //scaling factor L.columnSpacing is multiplied by length of RoI
            //scaling factor pinWheel count is divided by area of counting
            pinDensity = pinCount*L.columnSpacing*L.columnSpacing/roIarea;
            pinData << " " << pinCount << " " << pinDensity <<  std::endl;
            cout<<"pattern frequency " << L.patternFrequency << "  columnSpacing " << L.columnSpacing << " pinWheel density " << pinDensity << std::endl;
            cout << "max arg of psi  " << L.maxVal(psiphase) << " min arg of psi " << L.minVal(psiphase) <<  " iteration " << i <<endl;
            cerr << "max arg of psi  " << L.maxVal(psiphase) << " min arg of psi " << L.minVal(psiphase) <<  " iteration " << i <<endl;
            maxData << psiMax << " " << psiMin  << std::endl;
            cout << "max val of abs(psi)  " << psiMax << " min val of abs(psi) " << psiMin <<  " iteration " << std::endl;
            cerr << "max val of abs(psi)  " << psiMax << " min val of abs(psi) " << psiMin <<  " iteration " << std::endl;
#ifdef COMPILE_PLOTTING
            array<float,3> cl_a = {0.0f, 1.0f, 0.0f};
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
                    (*pW)[0] -= hexWidth*0.60;
                    (*pW)[1] += hexWidth*0.60;
                    vtx += morph::Vector<float, 3>({1,0,0});
                    //v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW, vtx, sz, 0.01f, cl_c, 6));
                    (*pW)[0] += hexWidth+0.10f;
                    (*pW)[1] += hexWidth*0.0;
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW, vtx, sz/2.0f, 0.01f, cl_b, 6));
                }
                std::cout << "just after pW in  lShowPinCentres" << std::endl;
                for (pW1 = pWCoords1.begin(); pW1 != pWCoords1.end(); pW1++) {
                    morph::Vector<float,3> vtx = *pW1;
                    (*pW1)[0] -= hexWidth*0.60;
                    (*pW1)[1] += hexWidth*0.60;
                    vtx += morph::Vector<float, 3>({1,0,0});
                    //v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW1, vtx, sz, 0.01f, cl_c, 6));
                    (*pW1)[0] += hexWidth+0.10f;
                    (*pW1)[1] += hexWidth*0.0;
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW1, vtx, sz/2.0f, 0.01f, cl_a, 6));
                }
                std::cout << "just after pWi1 in  lShowPinCentres" << std::endl;
                for (pW3 = pWCoords3.begin(); pW3 != pWCoords3.end(); pW3++) {
                    morph::Vector<float,3> vtx = *pW3;
                    (*pW3)[0] -= hexWidth*0.60;
                    (*pW3)[1] += hexWidth*0.60;
                    vtx += morph::Vector<float, 3>({1,0,0});
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW3, vtx, sz, 0.01f, cl_c, 6));
                    (*pW3)[0] += hexWidth+0.10f;
                    (*pW3)[1] += hexWidth*0.0;
                    v1.addVisualModel (new morph::PolygonVisual (v1.shaderprog, offset2, *pW3, vtx, sz, 0.01f, cl_c, 6));
                    pwCount ++;
                }


                std::cout << "pinWheel count via pWCoords3 " << pwCount << std::endl;
            }
            hgv1->updateData (&psiphase);
            hgv1->clearAutoscaleColour();

            hgv2->updateData (&psir);
            hgv2->clearAutoscaleColour();

            hgv3->updateData (&phiphase);
            hgv3->clearAutoscaleColour();


            graphX = L.xs;
            graphY = L.ys;
            /*
            for (int i=0;i<nbins/2;i++){
                std::cout << " frequency " << graphX[i] << " power " << graphY[i] << std::endl;
            }
            */
            int nsamp = nbins-10;
            float xmax = nsamp;

            arma::vec xfit(nsamp);
            graphX2.resize(nsamp,0);
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

            graphX3[0] = L.patternFrequency;
            graphX3[1] = L.patternFrequency;

            gvPinDensity->update (graphX, graphY, 0);
            gvPinDensity->update (graphX2, graphY2, 1);
            gvPinDensity->update (graphX3, graphY3, 2);

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

#endif
        } //end of if on numprint


     } //end of numsteps loop
//cout << " just after time step i = " << i << endl;

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
    outdata.add_contained_vals("psiR",psir); //Laplacian of psi
    outdata.add_contained_vals("psiPhase",psiphase); //main variable
    outdata.add_contained_vals("phiR",phir); //Laplacian of psi
    outdata.add_contained_vals("phiPhase",phiphase); //main variable
    //outdata.add_contained_vals("psiImg", psiImg); //psi as a matrix
    //data.add_val ("/g", g); //g parameter
    //data.add_val ("/epsilon", epsilon); //epsilon paramater
    cout << " just after writing data "  << endl;
//
/*
#ifdef COMPILE_PLOTTING
    cout << "Ctrl-c or press x in graphics window to exit.\n";
    v1.keepOpen();
#endif
*/
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
