#include <morph/Visual.h>
#include <morph/VisualDataModel.h>
#include <morph/HexGridVisual.h>
#include <morph/HexGrid.h>
#include <morph/ReadCurves.h>
#include <morph/tools.h>
#include <utility>
#include <iostream>
#include <fstream>
#include <cmath>
#include <morph/Scale.h>
#include <morph/Vector.h>

using namespace std;
using morph::Visual;
using morph::VisualDataModel;
using morph::HexGrid;
using morph::HexGridVisual;
using morph::Tools;
using morph::HexDomainShape;
using morph::ReadCurves;
using morph::Scale;
using morph::Vector;

int main()
{
    int rtn = -1;

    Visual v(1000,1000,"Test window");
    v.zNear = 0.001;

    try {
        string pwd = Tools::getPwd();
        string curvepath = "./tests/trial.svg";
        if (pwd.substr(pwd.length()-11) == "build/tests") {
            curvepath = "./../tests/trial.svg";
        }
        ReadCurves r(curvepath);

        HexGrid hg(0.01, 3, 0, HexDomainShape::Boundary);
        hg.setBoundary (r.getCorticalPath());

        cout << hg.extent() << endl;

        cout << "Number of hexes in grid:" << hg.num() << endl;
        cout << "Last vector index:" << hg.lastVectorIndex() << endl;

        if (hg.num() != 1604) {
            rtn = -1;
        }

        vector<float> data;
        unsigned int nhex = hg.num();
        data.resize(nhex, 0.0);

        // Make some dummy data (a sine wave)
        for (unsigned int hi=0; hi<nhex; ++hi) {
            data[hi] = 0.5 + 0.5*std::sin(10*hg.d_x[hi]); // Range 0->1
        }
        cout << "Created " << data.size() << " floats in data" << endl;

        Vector<float, 3> offset = { 0.0, 0.0, 0.0 };
        unsigned int gridId = v.addVisualModel (new HexGridVisual<float>(v.shaderprog, &hg, offset, &data));
        cout << "Added HexGridVisual with gridId " << gridId << endl;

        // Divide existing scale by 1:
        float newGrad = static_cast<VisualDataModel<float>*>(v.getVisualModel(gridId))->zScale.getParams(0)/1.0;
        // Set this in a new zscale object:
        Scale<float> zscale;
        zscale.setParams (newGrad, 0);
        // And set it back into the visual model:
        static_cast<VisualDataModel<float>*>(v.getVisualModel(gridId))->setZScale (zscale);


        vector<float> data1;
        data1.resize(nhex, 0.0);

        // Make some dummy data (a sine wave)
        for (unsigned int hi=0; hi<nhex; ++hi) {
            data1[hi] = 0.5 + 0.5*std::sin(10*hg.d_y[hi]); // Range 0->1
        }
        cout << "Created " << data1.size() << " floats in data" << endl;

        offset = { 1.0, 0.0, 1.0 };

        unsigned int gridId1 = v.addVisualModel (new HexGridVisual<float>(v.shaderprog, &hg, offset, &data1));
        cout << "Added another HexGridVisual with gridId " << gridId1 << endl;


        // Make some dummy data (a sine wave)
        for (unsigned int hi=0; hi<nhex; ++hi) {
            data1[hi] = 0.5 + 0.5*std::sin(10*hg.d_y[hi]*hg.d_x[hi]); // Range 0->1
        }
        cout << "Created " << data1.size() << " floats in data" << endl;

        offset = { 0.0, -0.8, 0.5 };

        unsigned int gridId2 = v.addVisualModel (new HexGridVisual<float>(v.shaderprog, &hg, offset, &data1));
        cout << "Added another HexGridVisual with gridId " << gridId2 << endl;

        v.render();

        while (v.readyToFinish == false) {
            glfwWaitEventsTimeout (0.018);
            v.render();
        }

    } catch (const exception& e) {
        cerr << "Caught exception reading trial.svg: " << e.what() << endl;
        cerr << "Current working directory: " << Tools::getPwd() << endl;
        rtn = -1;
    }


    return rtn;
}
