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
        float xspan = conf.getFloat("xspan",5.0);
#else
        double xspan = conf.getDouble("xspan",5.0);
#endif
    int scale = conf.getInt("scale",8);
    int fov = conf.getInt("fov",50);
    bool LfixedSeed = conf.getBool("LfixedSeed",0);
    //bool overwrite_logs = conf.getBool("overwrite_logs",true);
    bool skipMorph  = conf.getBool("skipMorph",false);
    int  iPolygon = conf.getInt("iPolygon", 2);
    unsigned int numpoints = conf.getInt("numpoints",41);
    vtxVisual* cv;
    int framenum = 0;

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
//end of integration after solving on polygonal regions via ksSolver

    /*
     * now create an array of morphed regions, first morph
     */

// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(0);
    cout << "just after setting curved boundaries morph 1 " << M.curvedBoundary.size()<<endl;
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
    // to draw the tessellation
// first convert M.vCoords into a morph::Vect
    cout << "just beofere creating vtxVector" << endl;
    vtxVector.resize(M.vCoords.size());
    cout << "just after creating vtxVector" << endl;
    for (unsigned int j=0;j<numpoints;j++) {
        vtxVector[j].resize(0);
        for (unsigned idx = 0; idx<M.vCoords[j].size(); idx++) {
            vtxVector[j].push_back(M.hGeo->point2vect3(M.vCoords[j][idx]));
        }
    }
    cout << "after filling vtxVector " << endl;
    // now instantiate vtxVisual

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
    v2->render();
//begin second morphing
// now set the boundaries for the regions, stored in curvedBoundary
    M.populateBoundCurve(1);
    cout << "just after setting curved boundaries morph 2 " << M.curvedBoundary.size()<<endl;
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
        // to draw the tessellation
// first convert M.vCoords into a morph::Vect
    cout << "just beofere creating vtxVector" << endl;
    vtxVector.resize(M.vCoords.size());
    cout << "just after creating vtxVector" << endl;
    for (unsigned int j=0;j<numpoints;j++) {
        vtxVector[j].resize(0);
        for (unsigned idx = 0; idx<M.vCoords[j].size(); idx++) {
            vtxVector[j].push_back(M.hGeo->point2vect3(M.vCoords[j][idx]));
        }
    }
    cout << "after filling vtxVector " << endl;
    // now instantiate vtxVisual
    cv = new vtxVisual(v3->shaderprog, v3->tshaderprog, vtxVector, M.ds, M.ds);
    cout << " after creating vtxVisual " << endl;
    cv->finalize();
    cout << " after vtxVisual finalize " << endl;
    cv->addLabel("Tessellation morph 2", {0.4f, 1.4f, 0.0f});
    v3->addVisualModel(cv);
    v3->render();
    frame << "log/agent/";
    frame.width(4);
    frame.fill('0');
    frame << framenum++;
    frame << ".png";
    savePngs (logpath, "tessellation2", 0, *v3);
    v3->setCurrent();
    v3->render();
    return 0;
};
