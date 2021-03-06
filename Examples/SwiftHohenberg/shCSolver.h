/*
  shSolver class
 * Author: John Brooke
 *
 * Date 2019/10
 *
 * creates a hexGrid given a boundary curve and solves
 * the KS equations
 *
 */
#include <morph/tools.h>
#include <morph/ReadCurves.h>
#include <morph/HdfData.h>
#include <morph/BezCurve.h>
#include <morph/BezCurvePath.h>
#include <morph/BezCoord.h>
#include <morph/HexGrid.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <iomanip>
#include <math.h>
#include <random>
#include <algorithm>
#include <hdf5.h>
#include <unistd.h>
#include <bits/stdc++.h>
#include <iostream>
#include <sys/stat.h>
#include <sys/types.h>
#ifndef PI
#define PI 3.1415926535897932f
#endif
using morph::HexGrid;
using morph::HdfData;
using morph::ReadCurves;
using morph::Tools;
using morph::BezCurve;
using morph::BezCurvePath;
using morph::BezCoord;
using namespace std;

#ifndef NO_N
#define NO_N
#define nE(hi) (this->Hgrid->d_ne[hi])
#define HAS_nE(hi) (this->Hgrid->d_ne[hi] == -1 ? false : true)

#define nW(hi) (this->Hgrid->d_nw[hi])
#define HAS_nW(hi) (this->Hgrid->d_nw[hi] == -1 ? false : true)

#define nNE(hi) (this->Hgrid->d_nne[hi])
#define HAS_nNE(hi) (this->Hgrid->d_nne[hi] == -1 ? false : true)

#define nNW(hi) (this->Hgrid->d_nnw[hi])
#define HAS_nNW(hi) (this->Hgrid->d_nnw[hi] == -1 ? false : true)

#define nSE(hi) (this->Hgrid->d_nse[hi])
#define HAS_nSE(hi) (this->Hgrid->d_nse[hi] == -1 ? false : true)

#define nSW(hi) (this->Hgrid->d_nsw[hi])
#define HAS_nSW(hi) (this->Hgrid->d_nsw[hi] == -1 ? false : true)
#endif

class shSolver
{
public:
// list of objects visible to member functions
    int scale; //sets hext to hex distance
    FLT xspan; //width of the hexGrid
    int n; //number of hexes in the HexGrid
    FLT ds; //distance used in numerical approximations
    FLT nnInitialOffset = 1.0;
    FLT ccInitialOffset = 2.5;
    FLT boundaryFalloffDist = 0.024;
    FLT sigma; //for the Gaussian convolution
    FLT gNorm; //denominator for convolution integrals
    pair<FLT,FLT> seedPoint;
    BezCurvePath<float> bound;
    string logpath;
    vector<vector<int> > N; // hex neighbourhood
    vector<complex<FLT>> psi, phi; //hold the field values for each hex
    vector<FLT> kerneldata;
    morph::HexGrid* Hgrid;
    morph::HexGrid* kernel;
// empty constructor
    shSolver(){};
// constructor with HexGrid passed in
    shSolver (morph::HexGrid*  Hgrid, std::string logpath) {
        this->Hgrid = Hgrid;
        this->seedPoint = this->Hgrid->computeCentroid(this->Hgrid->hexen);
        this->logpath = logpath;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        this->seedPoint = this->Hgrid->computeCentroid(this->Hgrid->hexen);
        this->ds = this->Hgrid->hexen.begin()->d; //hex-hex distance
        n = 0;
        n = this->Hgrid->num();
        afile << " max x " << this->Hgrid->getXmax(0.0) << " min x " << this->Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << this->Hgrid->num() << endl;
        afile << "after creating HexGri ds = "<< this->ds << endl;
        cout << "after creating HexGri ds = "<< this->ds << endl;
        afile << "seed point.x " << seedPoint.first << " seed point.y " << seedPoint.second << endl;
        n = this->Hgrid->num();
        afile << " max x " << this->Hgrid->getXmax(0.0) << " min x " << this->Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in boundary constructor" << " n = " << n <<endl;
      // check the order numbering in hexen
        this->N.resize(n);
        this->setHexType();
        this->psi.resize(n);
        this->phi.resize(n);
        cout << " end of shSolver from file " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
    }; // end of shSolver constructor

// constructor with radius passed in for solving on radial boundaries
    shSolver (int scale, FLT xspan, string logpath, float radius, pair<float, float> seedPoint) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << Hgrid->num() << endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        afile << "after creating HexGrid with radius "<< radius << endl;
        Hgrid->setCircularBoundary(radius, seedPoint, false);
        afile << "after setting boundary on  H " << Hgrid->num() << " centroid.x " << Hgrid->boundaryCentroid.first << " centroid.y " << Hgrid->boundaryCentroid.second << endl;
        afile << "after setting boundary on  H " << Hgrid->num() << endl;
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in circular constructor" << " n = " << n <<endl;
        N.resize(n);
        this->setHexType();
        this->psi.resize(n);
        this->phi.resize(n);
        cout << " end of shSolver circle radius " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
    }; // end of shSolver constructor

// Constructor with boundary passed in
    shSolver (int scale, FLT xspan, string logpath, BezCurvePath<float> bound, pair<FLT,FLT> seedPoint) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << Hgrid->num() << endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        Hgrid->setBoundary(bound,false);
        reverse_y();
        afile << "after setting boundary on  H " << Hgrid->num() << " centroid.x " << Hgrid->boundaryCentroid.first << " centroid.y " << Hgrid->boundaryCentroid.second << endl;
        afile << "seed point.x " << seedPoint.first << " seed point.y " << seedPoint.second << endl;
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in boundary constructor" << " n = " << n <<endl;
      // check the order numbering in hexen
        N.resize(n);
        this->setHexType();
        this->psi.resize(n);
        this->phi.resize(n);
        cout << " end of shSolver bezCurvePath " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
    }; // end of shSolver constructor



// Constructor for parallelogram domain
    shSolver (int scale, FLT xspan, string logpath) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        FLT pspan = xspan/3.0;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Parallelogram);

        int rextent =  floor(0.5 * pspan / this->ds) + 1;
        int gextent = rextent;
        cout << "before setting parallelogram boundary " << endl;
        Hgrid->setParallelogramBoundary(rextent, gextent);
        cout << "after setting parallelogram boundary " << endl;
        this->n = Hgrid->num();
        std::cout << "length of rows " << Hgrid->d_rowlen << " number of rows " << Hgrid->d_numrows << std::endl;
        std::cout << " size of HexGrid " << Hgrid->d_size << std::endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        N.resize(n);
        afile << " after N.resize()" << endl;
        this->setHexType();
        this->psi.resize(n);
        afile << "after psi.resize()" << endl;
        this->phi.resize(n);
        afile << "after psi.resize()" << endl;
        //now set up kernel
        this->sigma = this->xspan / 37.5;
        this->gNorm = 1.0 /(2.0*PI*sigma*sigma);
        this->kernel = new morph::HexGrid(this->ds, 20.0*sigma, 0, morph::HexDomainShape::Boundary);
        this->kernel->setCircularBoundary(6.0*sigma);
        this->kerneldata.resize(this->kernel->num(),0.0);
        cout << " end of shSolver parallelogram " << endl;

    }; // end of shSolver constructor


// Constructor for rectangular domain
    shSolver (int scale, FLT xspan, string logpath, float x, float y) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        n = 0;
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Rectangle);

        n = Hgrid->num();
        afile << "after setting rectangle boundary " << endl;
        Hgrid->setRectangularBoundary(x, y);
        afile << "after setting rectangle boundary " << endl;
        Hgrid->populate_d_neighbours();
        afile << "after populating d_neighbours " << endl;
        afile << "before filling H " << n << endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        //reverse_y();
      // check the order numbering in hexen
        N.resize(n);
        afile << " after N.resize()" << endl;
        //this->setHexType();
        this->psi.resize(n);
        afile << "after psi.resize()" << endl;
        this->phi.resize(n);
        afile << "after psi.resize()" << endl;
        cout << " end of shSolver bezCurvePath " << endl;

    }; // end of shSolver constructor

    //method to set the type of hexes
    void setHexType() {
        cout << "in setHexType" << endl;
        for (auto &h : this->Hgrid->hexen){
            this->N[h.vi].resize(6);
            if (!HAS_nE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][0] = h.vi;
            }
            else {
                this->N[h.vi][0] = Hgrid->d_ne[h.vi];
            }

            if (!HAS_nNE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][1] = h.vi;
            }
            else {
                this->N[h.vi][1] = Hgrid->d_nne[h.vi];
            }

            if (!HAS_nNW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][2] = h.vi;
            }
            else {
                this->N[h.vi][2] = Hgrid->d_nnw[h.vi];
            }

            if (!HAS_nW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][3] = h.vi;
            }
            else {
                this->N[h.vi][3] = Hgrid->d_nw[h.vi];
            }

            if (!HAS_nSW(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][4] = h.vi;
            }
            else {
                this->N[h.vi][4] = Hgrid->d_nsw[h.vi];
            }

            if (!HAS_nSE(h.vi)) {
                h.setBoundaryHex();
                this->N[h.vi][5] = h.vi;
            }
            else {
                this->N[h.vi][5] = Hgrid->d_nse[h.vi];
            }
        } //end of loop over HexGri this->setHexType();
    } // end of setHexType

    void setKerneldata() {
    // Once-only parts of the calculation of the Gaussian.
        FLT one_over_sigma_root_2_pi = 1 / sigma * 2.506628275;
        FLT two_sigma_sq = 2.0f * sigma * sigma;
        // Gaussian dist. result, and a running sum of the results:
        FLT gauss = 0;
        FLT sum = 0;
        for (auto& k : this->kernel->hexen) {
            // Gaussian profile based on the hex's distance from centre, which is
            // already computed in each Hex as Hex::r
            gauss = (one_over_sigma_root_2_pi * std::exp ( -(k.r*k.r) / two_sigma_sq ));
            kerneldata[k.vi] = gauss;
            sum += gauss;
        }
        // Renormalise
        for (auto& k : this->kernel->hexen) { this->kerneldata[k.vi] /= sum; }
    }

    //method to produce a Gaussian kernel
    vector<FLT> gaussConvolve(vector<FLT> invector) {
        vector <FLT> convolved;
        convolved.resize(invector.size(), 0.0);

        //now convolve
        this->Hgrid->convolve(*(this->kernel), this->kerneldata, invector, convolved);
        return convolved;
    }

// method to calculate the Laplacian
    vector<complex<FLT>> getLaplacian(vector<complex<FLT>> Q, FLT dx) {
        FLT overds = 1./(1.5*29.0*29.0*dx*dx);
        vector<complex<FLT>> L(n,0.);
        for(auto &h : this->Hgrid->hexen){
         int i = int(h.vi);
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.0f*Q[i])*overds;
        }
        return L;
    }

    vector<complex<FLT>> chemoTaxis(vector<complex<FLT>> Q, vector<complex<FLT>> P, FLT dx) {
        vector<complex<FLT>> cT(n,0.);
        FLT overds = 1.0f/(1.5f*29.0f*29.0f*dx*dx);
        FLT denom = 1.0f/2.0f;
        for (auto &h : Hgrid->hexen) {
            unsigned int i = h.vi;
        // finite volume method Lee et al. https://doi.org/10.1080/00207160.2013.864392
            complex<FLT> dr0Q = (Q[N[i][0]]+Q[i]) * denom;
            complex<FLT> dg0Q = (Q[N[i][1]]+Q[i]) * denom;
            complex<FLT> db0Q = (Q[N[i][2]]+Q[i]) * denom;
            complex<FLT> dr1Q = (Q[N[i][3]]+Q[i]) * denom;
            complex<FLT> dg1Q = (Q[N[i][4]]+Q[i]) * denom;
            complex<FLT> db1Q = (Q[N[i][5]]+Q[i]) * denom;

            complex<FLT> dr0P = P[N[i][0]]-P[i];
            complex<FLT> dg0P = P[N[i][1]]-P[i];
            complex<FLT> db0P = P[N[i][2]]-P[i];
            complex<FLT> dr1P = P[N[i][3]]-P[i];
            complex<FLT> dg1P = P[N[i][4]]-P[i];
            complex<FLT> db1P = P[N[i][5]]-P[i];


            cT[i] = (dr0Q*dr0P+dg0Q*dg0P+db0Q*db0P+dr1Q*dr1P+dg1Q*dg1P+db1Q*db1P)*overds;

        } //matches for on i
        return cT;
    } //end of function chemoTaxis

  // function to compute the derivative of the Swift Hohenberg equation
     void compute_dpsidt(vector<complex<FLT>>& inPsi, vector<complex<FLT>>& dpsidt, vector<complex<FLT>> inPhi, FLT epsilon, FLT g, bool lRange=false) {
        vector<complex<FLT>> lapPhi(this->n,0);
        vector <FLT> psimodsq;
        int size = inPsi.size();
        //std::cerr << "in compute_dpi size  " << size <<  std::endl;
        psimodsq.resize(size, 0.0);
        FLT zmod = 0;
        for (int i=0; i<size; i++) {
            zmod = abs(inPsi[i]);
            psimodsq[i] = zmod*zmod;
        }
        //std::cerr << "just before gaussConvolve " << std::endl;
        //vector <FLT> integral = this->gaussConvolve(psimodsq);
        //std::cerr << "just after gaussConvolve " << std::endl;
        lapPhi = getLaplacian(inPhi,this->ds);

        if (lRange == true) {
            for (int h=0; h < this->n; h++) {
                dpsidt[h] = -2.0f*this->phi[h] - lapPhi[h] + (epsilon - 1.0f)*inPsi[h] +  (1.0f - g) * psimodsq[h]*inPsi[h];
          //      dpsidt[h] += (2.0f-g) * integral[h] * inPsi[h] / this->gNorm;
            }
        }
        else {
            for (int h=0; h < this->n; h++) {
                dpsidt[h] = -2.0f*this->phi[h] - lapPhi[h] + (epsilon - 1.0f)*inPsi[h] +  (1.0f - g) * psimodsq[h]*inPsi[h];
            }
        }
    }//end of method compute_dpsidt


     void compute_dpsidtH(vector<complex<FLT>>& inPsi, vector<complex<FLT>>& dpsidt, vector<complex<FLT>> inPhi, FLT epsilon, FLT g) {
        vector<complex<FLT>> lapPsi(this->n,0);
        lapPsi = getLaplacian(inPsi,this->ds);
        for (int di=0; di < this->n; di++) {
            dpsidt[di] =  lapPsi[di];
        }
    }

    void compute_phi (vector<complex<FLT>>& inPsi, vector<complex<FLT>>& oldPsi) {
        vector<complex<FLT>> newLaplacian;
        vector<complex<FLT>> oldLaplacian;
        newLaplacian = getLaplacian(inPsi,this->ds);
        oldLaplacian = getLaplacian(oldPsi,this->ds);
        for (unsigned int i=0; i<inPsi.size(); i++) {
            this->phi[i] = 0.5f*(oldLaplacian[i] + newLaplacian[i]);
        }
    }


  //function to time step Neumann b.c.s for the Swift-Hohenberg equation
    void step(FLT dt, FLT epsilon, FLT g, vector<complex<FLT>> oldPsi)
    {
        /* Set up boundary conditions with ghost points
        //cout << " in time step before ghost points" << endl;
        for(auto &h : Hgrid->hexen)
        {
            if (h.boundaryHex())
            {
                for(int j=0;j<6;j++)
                {
                    int i = int(h.vi);
                    if(N[h.vi][j] == i) {
                        this->psi[N[h.vi][j]] = this->psi[h.vi];
                        this->phi[N[h.vi][j]] = this->phi[h.vi];
                    }
                }
            }
        }
        //cerr << "in step after setting b.c.s " << std::endl;
        */
        // 2. Do integration of psi
        {
            // Runge-Kutta integration for psi. This time, I'm taking
            // ownership of this code and properly understanding it.

            // Ntst: "A at a test point". Ntst is a temporary estimate for A.
            vector<complex<FLT>> Ntst(this->n, 0.0);
            vector<complex<FLT>> dpsidt(this->n, 0.0);
            vector<complex<FLT>> K1(this->n, 0.0);
            vector<complex<FLT>> K2(this->n, 0.0);
            vector<complex<FLT>> K3(this->n, 0.0);
            vector<complex<FLT>> K4(this->n, 0.0);

            /*
             * Stage 1
             */
            this->compute_dpsidt (this->psi, dpsidt, this->phi, epsilon, g);
            for (int h=0; h< this->n; ++h) {
                K1[h] = dpsidt[h] * dt;
                Ntst[h] = this->psi[h] + K1[h] * 0.5f ;
            }

            /*
             * Stage 2
             */
            this->compute_dpsidt (Ntst, dpsidt, this->phi, epsilon, g);
            for (int h=0; h< this->n; ++h) {
                K2[h] = dpsidt[h] * dt;
                Ntst[h] = this->psi[h] + K2[h] * 0.5f;
            }

            /*
             * Stage 3
             */
            this->compute_dpsidt (Ntst, dpsidt, this->phi, epsilon, g);
            for (int h=0; h < this->n; ++h) {
                K3[h] = dpsidt[h] * dt;
                Ntst[h] = this->psi[h] + K3[h];
            }

            /*
             * Stage 4
             */
            this->compute_dpsidt (Ntst, dpsidt, this->phi, epsilon, g);
            for (int h=0; h < this->n; ++h) {
                K4[h] = dpsidt[h] * dt;
            }


            /*
             * Final sum together. This could be incorporated in the
             * for loop for Stage 4, but I've separated it out for
             * pedagogy.
             */
            for (int h=0;h<this->n;h++) {
                this->psi[h] += ((K1[h] + 2.0f * (K2[h] + K3[h]) + K4[h])/(FLT)6.0);
            }
        }

    // 3. Do integration of phi
       compute_phi(this->psi, oldPsi);
    }//end step

  //function to time step periodic b.c.s for the heat equation
    void stepHeat(FLT dt, FLT epsilon, FLT g, vector<complex<FLT>> oldPsi)
    {
        // Set up boundary conditions with ghost points
        //cout << " in time step before ghost points" << endl;
        for(auto &h : Hgrid->hexen)
        {
            if (h.boundaryHex())
            {
                for(int j=0;j<6;j++)
                {
                    int i = int(h.vi);
                    if(N[h.vi][j] == i) {
                        this->psi[N[h.vi][j]] = this->psi[h.vi];
                        this->phi[N[h.vi][j]] = this->phi[h.vi];
                    }
                }
            }
        }
        // 1. Do integration of psi
        // Runge-Kutta integration for psi. This time, I'm taking
        // ownership of this code and properly understanding it.

        // Ntst: "A at a test point". Ntst is a temporary estimate for A.
        vector<complex<FLT>> Ntst(this->n, 0.0);
        vector<complex<FLT>> dpsidt(this->n, 0.0);
        vector<complex<FLT>> K1(this->n, 0.0);
        vector<complex<FLT>> K2(this->n, 0.0);
        vector<complex<FLT>> K3(this->n, 0.0);
        vector<complex<FLT>> K4(this->n, 0.0);

        /*
            * Stage 1
        */
        this->compute_dpsidtH (this->psi, dpsidt, this->phi, epsilon, g);
        for (int di=0; di< this->n; ++di) {
            K1[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K1[di] * 0.5f ;
        }

        /*
         * Stage 2
        */
        this->compute_dpsidtH (Ntst, dpsidt, this->phi, epsilon, g);
        for (int di=0; di< this->n; ++di) {
            K2[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K2[di] * 0.5f;
        }

        /*
         * Stage 3
         */
        this->compute_dpsidtH (Ntst, dpsidt, this->phi, epsilon, g);
        for (int di=0; di < this->n; ++di) {
            K3[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K3[di];
        }

        /*
         * Stage 4
         */
        this->compute_dpsidtH (Ntst, dpsidt, this->phi, epsilon, g);
        for (int di=0; di < this->n; ++di) {
            K4[di] = dpsidt[di] * dt;
        }


        /*
         * Final sum together. This could be incorporated in the
         * for loop for Stage 4, but I've separated it out for
         * pedagogy.
         */
        for (int di=0;di<this->n;di++) {
            this->psi[di] += ((K1[di] + 2.0f * (K2[di] + K3[di]) + K4[di])/(FLT)6.0);
        }

        // 2. Do integration of phi
        compute_phi(this->psi, oldPsi);
    }//end step

    void reverse_y ()
    {
        for (auto &h : this->Hgrid->hexen)
        {
            FLT temp = h.y;
            if (temp != 0.0)
            {
                h.y = -temp;
                this->Hgrid->d_y[h.vi] = -temp;
            }
        }
    }


    vector<FLT> complexZero(vector<FLT> invector) {
        vector<FLT> result;
        if (static_cast<int> (invector.size()) != this->n) {
            //cerr << " compleZero invector " << invector.size() << " not equal to this->n " << this->n << endl;
            std::exit(0);
        }
        int count = 0;
        for(auto h: this->Hgrid->hexen) {
            bool test = true;
            if (h.boundaryHex()) continue;
            for (int i=0; i<6; i++) {
                if (invector[h.vi] > invector[N[h.vi][i]] && N[h.vi][i] > 0) {
                    test = false;
                    break;
                }
            }
            if (test == true) {
                result.push_back(invector[h.vi]);
                count++;
            }
        }
        cout << " number of pinwheels " << count << " out of " << this->n << " hexes " << endl;
        return result;
    }

    vector<bool> cZero(vector<FLT> invector, FLT epsilon){
        vector<bool> result;
        result.resize(invector.size(), false);
        int count = 0;
        for (auto h: this->Hgrid->hexen) {
            if (invector[h.vi] < epsilon) {
                result[h.vi] = true;
                count++;
            }
        }
        cout << " number of pinwheels " << count << " out of " << this->n << " hexes " << endl;
        return result;
    }





    vector<bool> contour(vector<FLT> invector) {
        vector<bool> result;
        if (static_cast<int> (invector.size()) != this->n) {
            //cerr << " compleZero invector " << invector.size() << " not equal to this->n " << this->n << endl;
            std::exit(0);
        }
        result.resize(invector.size(), true);
        int count = 0;
        for(auto h: this->Hgrid->hexen) {
            for (int i=0; i<6; i++) {
                if (invector[h.vi] > invector[N[h.vi][i]] ) {
                    result[h.vi] = false;
                }
            }
            if (result[h.vi] == true) count++;
        }
        cout << " number of pinwheels " << count << " out of " << this->n << " hexes " << endl;
        return result;
    }
/* function to give r and theta relative to region centre
    pair <FLT,FLT> set_kS_polars(pair<FLT,FLT> centre){
        pair <FLT, FLT> result;
        result.first = 0.0;
        result.second = 0.0;
        FLT xav=0;
        FLT yav = 0;
        int hexcount = 0;
        FLT maxPhi = -10.0;
        FLT minPhi = 10.0;
        for (auto &h : this->Hgrid->hexen) {
            hexcount++;
            xav += h.x;
            yav += h.y;
        }
        if (hexcount != 0) {
            xav = xav / (hexcount*1.0);
            yav = yav / (hexcount*1.0);
        }
        else {
            cout << "in set ks_polars hexcount = 0" << endl;
        }
//go over the region and put the hexes into bins then average
        for (auto&  h : this->Hgrid->hexen) {
            FLT angle = 0;
            FLT dx = h.x;
            FLT dy = h.y;
            h.r = sqrt((dx - centre.first)*(dx - centre.first)
            + (dy - yav)*(dy - yav));
            if (dy >= centre.second) {
                angle =   atan2((dy - centre.second), (dx - centre.first));
                h.phi = angle;
            }
            else {
                angle =  2*PI + atan2((dy - centre.second), (dx - centre.first));
                h.phi = angle;
            }
            if (angle < minPhi) {
                minPhi = angle;
            }
            if (angle > maxPhi) {
                maxPhi = angle;
            }
        }
        cout << " set_kS_polars max phi " << maxPhi << " minPhi " << minPhi << endl;
        result.first = xav - centre.first ; // barycentre
        result.second = yav - centre.second;
		cout << "centre x "<< centre.first << " centre y " << centre.second << " centreMove.x " << result.first << " centreMove.y " << result.second <<endl;
        return result;
    } //end of function set_polars */

}; //end of class KSsolver
