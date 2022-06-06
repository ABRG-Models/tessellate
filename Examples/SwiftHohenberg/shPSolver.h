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
    FLT overds;
    FLT lengthScale;
    FLT nnInitialOffset = 1.0;
    FLT ccInitialOffset = 2.5;
    FLT boundaryFalloffDist = 0.024;
    FLT sigma; //for the Gaussian convolution
    FLT gNorm; //denominator for convolution integrals
    pair<FLT,FLT> seedPoint;
    BezCurvePath<float> bound;
    string logpath;
    vector<vector<int> > N; // hex neighbourhood
    vector<complex<FLT>> psi, phi; //hold the field values for each he
    vector<FLT> kernelRdata;
    vector<complex<FLT>> kernelCdata;
    vector<FLT> nonLocalR;
    vector<complex<FLT>> nonLocalC;
    morph::HexGrid* Hgrid;
    morph::HexGrid* kernel;
    vector<vector<int>> convIndex;
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
        this->n = 0;
        this->n = this->Hgrid->num();
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
        this->psi.resize(n);
        this->phi.resize(n);
        cout << " end of shSolver from file " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
    }; // end of shSolver constructor

// constructor with radius passed in for solving on radial boundaries
    shSolver (int scale, FLT xspan, string logpath, float radius, pair<float, float> seedPoint, FLT lengthScale) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        this->lengthScale = lengthScale;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        this->overds = 1.0/(this->lengthScale*this->lengthScale*ds*ds);
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
        this->setNoFlux();
        this->nonLocalR.resize(this->n, 0.0);
        this->nonLocalC.resize(this->n, 0.0);
        //now set up kernel
        this->sigma = 0.26; //based on freqency estimates of wavelength from simulations
        this->gNorm = 1.0 /(2.0*PI*sigma*sigma);
        this->kernel = new morph::HexGrid(this->ds, 10.0*sigma, 0, morph::HexDomainShape::Boundary);
        this->kernel->setCircularBoundary(3.0*sigma);
        this->kernelRdata.resize(this->kernel->num(),0.0);
        this->kernelCdata.resize(this->kernel->num());
        this->psi.resize(n);
        this->phi.resize(n);
        cout << " end of shSolver circle radius " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
    }; // end of shSolver constructor

// Constructor with boundary passed in
    shSolver (int scale, FLT xspan, string logpath, BezCurvePath<float> bound, pair<FLT,FLT> seedPoint, FLT lengthScale) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        this->lengthScale = lengthScale;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        this->overds = 1.0/(this->lengthScale*this->lengthScale*ds*ds);
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Boundary);
        n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "before filling H " << Hgrid->num() << endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        Hgrid->setBoundary(bound,false);
        reverse_y();
        afile << "after setting boundary on  H " << Hgrid->num() << " centroid.x " << Hgrid->boundaryCentroid.first << " centroid.y " << Hgrid->boundaryCentroid.second << endl;
        afile << "seed point.x " << seedPoint.first << " seed point.y " << seedPoint.second << endl;
        this->n = Hgrid->num();
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        afile << "after  filling H in boundary constructor" << " n = " << n <<endl;
      // check the order numbering in hexen
        N.resize(n);
        //this->setNoFlux();
        this->psi.resize(n);
        this->nonLocalR.resize(this->n, 0.0);
        this->nonLocalC.resize(this->n, 0.0);
        //now set up kernel
        this->sigma = this->xspan / 37.5;
        this->gNorm = 1.0 /(2.0*PI*sigma*sigma);
        this->kernel = new morph::HexGrid(this->ds, 10.0*sigma, 0, morph::HexDomainShape::Boundary);
        this->kernel->setCircularBoundary(3.0*sigma);
        this->kernelRdata.resize(this->kernel->num(),0.0);
        this->kernelCdata.resize(this->kernel->num());
        this->phi.resize(n);
        cout << " end of shSolver bezCurvePath " << " x seedPoint " << seedPoint.first << " y seedPoint " << seedPoint.second << endl;
    }; // end of shSolver constructor

// Constructor for parallelogram domain
    shSolver (int scale, FLT xspan, string logpath, FLT lengthScale) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        this->lengthScale = lengthScale;
        FLT pspan = xspan/3.0;
        ofstream afile (this->logpath + "/ksdebug.out",ios::app );
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        this->overds = 1.0/(this->lengthScale*this->lengthScale*ds*ds);
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Parallelogram);

        int rextent =  floor(0.5 * pspan / this->ds) + 1;
        int gextent = rextent;
        afile << "before setting parallelogram boundary " << endl;
        Hgrid->setParallelogramBoundary(rextent, gextent);
        afile << "after setting parallelogram boundary " << endl;
        this->n = Hgrid->num();
        cerr << "this->n " << this->n << " this->d_size " << Hgrid->d_size <<std::endl;
        std::cout << "length of rows " << Hgrid->d_rowlen << " number of rows " << Hgrid->d_numrows << std::endl;
        std::cout << " size of d_vexto " << Hgrid->d_size << " size of HexGrid " << Hgrid->num() << std::endl;
        afile << "after creating HexGrid ds =  " << this->ds << endl;
        afile << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << endl;
        this->psi.resize(this->n);
        afile << "after psi.resize()" << endl;
        this->phi.resize(this->n);
        this->nonLocalR.resize(this->n, 0.0);
        this->nonLocalC.resize(this->n, 0.0);
        //now set up kernel
        this->sigma = 0.142857;
        this->gNorm = 1.0 /(2.0*PI*sigma*sigma);
        this->kernel = new morph::HexGrid(this->ds, 6.0*sigma, 0, morph::HexDomainShape::Boundary);
        this->kernel->setCircularBoundary(2.0*sigma);
        this->kernelRdata.resize(this->kernel->num(),0.0);
        this->kernelCdata.resize(this->kernel->num());
        afile << "after psi.resize()" << endl;
        cout << " end of shSolver parallelogram " << endl;

    }; // end of shSolver constructor


// Constructor for rectangular domain
    shSolver (int scale, FLT xspan, string logpath, FLT x, FLT y, FLT lengthScale) {
        this->scale = scale;
        this->xspan = xspan;
        this->logpath = logpath;
        this->bound = bound;
        this->seedPoint = seedPoint;
        this->lengthScale = lengthScale;
        FLT s = pow(2.0, this->scale-1);
        this->ds = 1.0/s;
        this->overds = 1.0/(this->lengthScale*this->lengthScale*ds*ds);
        Hgrid = new HexGrid(this->ds, this->xspan, 0.0, morph::HexDomainShape::Rectangle);
        std::cout << "after setting rectangle boundary " << std::endl;
        Hgrid->setRectangularBoundary(x, y);
        std::cout << "after setting rectangle boundary " << std::endl;
        Hgrid->populate_d_neighbours();
        this->n = Hgrid->num();
        N.resize(n);
        std::cout << "after populating d_neighbours n " << this->n << std::endl;
        std::cout << "after creating HexGrid ds =  " << this->ds << std::endl;
        std::cout << " max x " << Hgrid->getXmax(0.0) << " min x " << Hgrid->getXmin(0.0) << std::endl;
        //this->setNoFlux();
        std::cout << "before psi.resize()" << std::endl;
        std::complex<FLT> const zero(0.0, 0.0);
        this->psi.resize(n, zero);
        this->phi.resize(n, zero);
        std::cout << "after psi.resize()" << std::endl;
        this->nonLocalR.resize(this->n, 0.0);
        this->nonLocalC.resize(this->n, 0.0);
        //now set up kernel
        this->sigma = x / 12.0;
        this->gNorm = 1.0 /(2.0*PI*sigma*sigma);
        this->kernel = new morph::HexGrid(this->ds, 6.0*sigma, 0, morph::HexDomainShape::Boundary);
        this->kernel->setCircularBoundary(2.0f*sigma);
        this->kernelRdata.resize(this->kernel->num(),0.0);
        this->kernelCdata.resize(this->kernel->num());
        std::cout << " end of shSolver rectangular " << std::endl;

    }; // end of shSolver constructor


    void setKernelRdata() {
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
            kernelRdata[k.vi] = gauss;
            sum += gauss;
        }
        // Renormalise
        for (auto& k : this->kernel->hexen) { this->kernelRdata[k.vi] /= sum; }
    }


    void setKernelCdata() {
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
            kernelCdata[k.vi].real(gauss);
            sum += gauss;
        }
        // Renormalise
        for (auto& k : this->kernel->hexen) {
            FLT av = kernelCdata[k.vi].real() / sum;
            this->kernelCdata[k.vi].real(av);
            this->kernelCdata[k.vi].imag(0.0);
        }
    }

    void setConvolutionIndex() {
        this->convIndex.resize(this->n);
        this->Hgrid->convolveIndex(*(this->kernel),this->convIndex);
    }

    //method to set the non-local terms
    void setNonLocalR() {
        vector<FLT> psimodsq;
        psimodsq.resize(this->n, 0.0);
        //now convolve
        //std::cout << "before Hgrid->convolve setNonLocalR" << std::endl;
        for (int i=0; i<this->n-10; i++) {
            FLT zmod = abs(this->psi[i]);
            psimodsq[i] = zmod*zmod;
            this->nonLocalR[i] = 0.0;
        }
#pragma omp parallel for
        for (int i=0;i<this->n;i++) {
            for (auto kh: kernel->hexen) {
                if (convIndex[i][kh.vi] < 0 || convIndex[i][kh.vi] > this->n) {
                    //std::cout << "after Hgrid->convolve i " << i << " kh.vi " << kh.vi << " convIndex " << convIndex[i][kh.vi] << std::endl;
                    this->nonLocalR[i] = 0;
                    continue;
                }
                this->nonLocalR[i] += this->kernelRdata[kh.vi]*psimodsq[this->convIndex[i][kh.vi]];
            }
        }
        //this->Hgrid->convolve(*(this->kernel), this->kernelCdata, psimodsq, this->nonLocalR);
        //std::cout << "after Hgrid->convolve " << std::endl;
    }

    //method to set the non-local terms
    void setNonLocalC() {
        vector<complex<FLT>> psisq;
        psisq.resize(this->n, 0.0);
        //now convolve
       // std::cout << "before Hgrid->convolve setNonLocalC" << std::endl;
        for (int i=0; i<this->n-10; i++) {
            psisq[i] = this->psi[i]*this->psi[i];
            this->nonLocalC[i].real(0.0);
            this->nonLocalC[i].imag(0.0);
        }
#pragma omp parallel for
        for (int i=0;i<this->n-10;i++) {
            for (auto kh: kernel->hexen) {
                if (convIndex[i][kh.vi] < 0 || convIndex[i][kh.vi] > this->n) {
                    //std::cout << "after Hgrid->convolve i " << i << " kh.vi " << kh.vi << " convIndex " << convIndex[i][kh.vi] << std::endl;
                    this->nonLocalC[i].real(0.0);
                    this->nonLocalC[i].imag(0.0);
                    continue;
                }
                this->nonLocalC[i] += this->kernelRdata[kh.vi]*psisq[this->convIndex[i][kh.vi]];
            }
        }
        //this->Hgrid->convolve(*(this->kernel), this->kernelCdata, psisq, this->nonLocalC);
        //std::cout << "after Hgrid->convolve " << std::endl;
    }

    vector<vector<FLT>> fieldInRoI(vector<FLT> A, int iRoI, FLT halfWidth, std::pair<FLT,FLT> centre) {
        vector<vector<FLT>> roIField;
        std::pair<FLT,FLT> botLeft; //coordinates of bottom left
       // int iRoI = 2*floor(halfWidth/this->ds) + 1.0f; //rowlength of square
        roIField.resize(iRoI);
        for (int i=0; i<iRoI; i++) {
            roIField[i].resize(iRoI,0);
        }
        botLeft.first = centre.first - halfWidth;
        botLeft.second = centre.second - halfWidth;
        int n = A.size();
        std::cout << "size of A " << n << " size of iRoI " << iRoI << std::endl;
        FLT maxV = -1e9;
        FLT minV = +1e9;
        for(int i=0; i<n; i++){
            if(maxV<fabs(A[i])){
                maxV = fabs(A[i]);
            } // THIS WAY ENSURE THAT ZERO DIFF ALWAYS MAPS TO VAL OF 0.5
            if(minV>A[i]){ minV = A[i]; }
        }
        FLT scale = 1./(2*maxV);

        //algortithm needs start of each row, first row begins bottom left
        std::list<morph::Hex>::iterator blh  = this->Hgrid->findHexNearest(botLeft);
        std::list<morph::Hex>::iterator h;
        for (int i=0; i<iRoI; i++) {
            h = blh;
            for(int j=0; j<iRoI; j++) {
                roIField[i][j] = A[(*h).di]*scale;
                h = h->ne;
                //std::cout << " h xval " << this->Hgrid->d_x[(*h).di] << " h yval " << this->Hgrid->d_y[(*h).di] << " h.di " << (*h).di << " field " << A[(*h).di] << " roI field " << roIField[i][j] << std::endl;
            }
            if (i%2 == 0)
                blh = blh->nne;
            else
                blh = blh->nnw;
            std::cout << "blh i " << i << " x " << this->Hgrid->d_x[(*blh).vi] << " y " << this->Hgrid->d_y[(*blh).vi] << std::endl;
        }
        return roIField;
    }



// method to calculate the Laplacian
    vector<complex<FLT>> getLaplacian(vector<complex<FLT>> Q) {
        vector<complex<FLT>> L(n,0.);
        for(int di=0; di<this->n; di++){
            L[di]=(Q[Hgrid->d_ne[di]]+Q[Hgrid->d_nne[di]]+Q[Hgrid->d_nnw[di]]+Q[Hgrid->d_nw[di]]+Q[Hgrid->d_nsw[di]]+Q[Hgrid->d_nse[di]]-6.0f*Q[di])*this->overds;
        }
        return L;
    }

    vector<complex<FLT>> chemoTaxis(vector<complex<FLT>> Q, vector<complex<FLT>> P) {
        vector<complex<FLT>> cT(n,0.);
        FLT denom = 1.0f/2.0f;
        for (int di=0; di<this->n; di++) {
        // finite volume method Lee et al. https://doi.org/10.1080/00207160.2013.864392
            complex<FLT> dr0Q = (Q[Hgrid->d_ne[di]]+Q[di]) * denom;
            complex<FLT> dg0Q = (Q[Hgrid->d_nne[di]]+Q[di]) * denom;
            complex<FLT> db0Q = (Q[Hgrid->d_nnw[di]]+Q[di]) * denom;
            complex<FLT> dr1Q = (Q[Hgrid->d_nw[di]]+Q[di]) * denom;
            complex<FLT> dg1Q = (Q[Hgrid->d_nsw[di]]+Q[di]) * denom;
            complex<FLT> db1Q = (Q[Hgrid->d_nse[di]]+Q[di]) * denom;

            complex<FLT> dr0P = P[Hgrid->d_ne[di]]-P[di];
            complex<FLT> dg0P = P[Hgrid->d_nne[di]]-P[di];
            complex<FLT> db0P = P[Hgrid->d_nnw[di]]-P[di];
            complex<FLT> dr1P = P[Hgrid->d_nw[di]]-P[di];
            complex<FLT> dg1P = P[Hgrid->d_nsw[di]]-P[di];
            complex<FLT> db1P = P[Hgrid->d_nse[di]]-P[di];


            cT[di] = (dr0Q*dr0P+dg0Q*dg0P+db0Q*db0P+dr1Q*dr1P+dg1Q*dg1P+db1Q*db1P)*this->overds;

        } //matches for on i
        return cT;
    } //end of function chemoTaxis

  /* function to compute the derivative
     void compute_dpsidt(vector<complex<FLT>>& inPsi, vector<complex<FLT>>& dpsidt, vector<complex<FLT>> inPhi, FLT epsilon, FLT g) {
        vector<complex<FLT>> lapPhi(this->n,0);
        lapPhi = getLaplacian(inPhi,this->ds);
        for (int di=0; di < this->n; di++) {
            FLT psimod = abs(phi[di]); //this should have been psi
            FLT psimodsq = psimod*psimod;
          dpsidt[di] = -2.0f*this->phi[di] - lapPhi[di] + (epsilon - 1.0f)*psi[di] +  (1.0f - g) * psimodsq*psi[di];
        }
    }
 */

  // function to compute the derivative of the Swift Hohenberg equation
     void compute_dpsidt(vector<complex<FLT>>& inPsi, vector<complex<FLT>>& dpsidt, FLT epsilon, FLT g, FLT k0) {
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
        lapPhi = getLaplacian(this->phi);

        for (int h=0; h < this->n; h++) {
            dpsidt[h] = -2.0f*k0*k0*this->phi[h] - lapPhi[h] + (epsilon - k0*k0*k0*k0)*inPsi[h] +  (1.0f - g) * psimodsq[h]*inPsi[h];
            dpsidt[h] -= (2.0f-g) * (this->nonLocalR[h] * inPsi[h] + 0.5f * this->nonLocalC[h] * std::conj(inPsi[h])) / this->gNorm;
        }
    }//end of method compute_dpsidt

    void compute_dpsidtH(vector<complex<FLT>>& inPsi, vector<complex<FLT>>& dpsidt,  FLT epsilon, FLT g) {
        vector<complex<FLT>> lapPsi(this->n,0);
        lapPsi = getLaplacian(inPsi);
        for (int di=0; di < this->n; di++) {
            dpsidt[di] =  lapPsi[di];
        }
    }

    void compute_phi (vector<complex<FLT>>& inPsi, vector<complex<FLT>>& oldPsi) {
        vector<complex<FLT>> newLaplacian;
        vector<complex<FLT>> oldLaplacian;
        newLaplacian = getLaplacian(inPsi);
        oldLaplacian = getLaplacian(oldPsi);
        for (unsigned int di=0; di<inPsi.size(); di++) {
            this->phi[di] = 0.5f*(oldLaplacian[di] + newLaplacian[di]);
        }
    }

    //method to set noFluxBC
    void setNoFlux() {
        cout << "in setNoFlux" << endl;
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
        } //end of loop over HexGri this->setNoFlux();
    } // end of setNoFlux


    //set periodic boundary conditions for a parallelogram grid
    void setPeriodic() {

        int di = 0;
        int offset = 0;
        //along bottom row first
        //bottom left corner
        Hgrid->d_nnw[0] = 2*Hgrid->d_rowlen-1;
        Hgrid->d_nw[0] = Hgrid->d_rowlen-1;
        Hgrid->d_nsw[0] = Hgrid->d_rowlen * (Hgrid->d_numrows-1);
        Hgrid->d_nse[0] = Hgrid->d_rowlen * (Hgrid->d_numrows-1) + 1;
        //bottom right hand corner
        Hgrid->d_nsw[Hgrid->d_rowlen-1] = Hgrid->d_rowlen*(Hgrid->d_numrows)-1;
        Hgrid->d_nse[Hgrid->d_rowlen-1] = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        Hgrid->d_ne[Hgrid->d_rowlen-1] = di;
        //rest of bottom row
        offset = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        for (unsigned int j=1; j<Hgrid->d_rowlen-1; j++) {
            Hgrid->d_nsw[j] = offset + j;
            Hgrid->d_nse[j] = offset + j+1;
        }
        std::cout << "bottom row  d_nsw[1] " << Hgrid->d_nsw[1] << std::endl;
        //top row now
        //top left can resuse offset
        Hgrid->d_nw[offset] = offset + Hgrid->d_rowlen - 1;
        Hgrid->d_nnw[offset] = Hgrid->d_rowlen-1;;
        Hgrid->d_nne[offset] = 0;
        //top right now
        offset = Hgrid->d_rowlen*(Hgrid->d_numrows) - 1;
        Hgrid->d_nnw[offset] = Hgrid->d_rowlen-2;
        Hgrid->d_nne[offset] = Hgrid->d_rowlen-1;
        Hgrid->d_ne[offset] = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        Hgrid->d_nse[offset] = Hgrid->d_rowlen*(Hgrid->d_numrows-2);
        //now the interior of the top row
        offset = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        for (unsigned int j=1; j<Hgrid->d_rowlen-1; j++) {
            di = offset + j;
            Hgrid->d_nnw[di] = j-1;
            Hgrid->d_nne[di] = j;
        }
        std::cout << "top row offset " << offset << " d_nnw " << Hgrid->d_nnw[offset] << std::endl;
        //now the b.c.s on the interior applied at start and end of each row
        for (unsigned int i=1; i<Hgrid->d_numrows-1; i++) {
            offset = i*Hgrid->d_rowlen;
            //left boundary
            Hgrid->d_nnw[offset] = offset + 2*Hgrid->d_rowlen - 1;
            Hgrid->d_nw[offset] = offset + Hgrid->d_rowlen - 1;
            //right boundary
            Hgrid->d_ne[offset+Hgrid->d_rowlen-1] = offset;
            Hgrid->d_nse[offset+Hgrid->d_rowlen-1] = offset - Hgrid->d_rowlen;
        }
        std::cout << "left col row 1 " << " d_nw " << Hgrid->d_nw[Hgrid->d_rowlen] << std::endl;
        std::cout << "right col row 1 " << " d_ne " << Hgrid->d_ne[2*Hgrid->d_rowlen-1] << std::endl;
    }

    //set boundary conditions for a rectangular domain with an the bottom row odd
    void setPeriodicOdd() {

        int di = 0;
        int offset = 0;
        //along bottom row first
        //bottom left corner
        Hgrid->d_nnw[0] = 2*Hgrid->d_rowlen-1;
        Hgrid->d_nw[0] = Hgrid->d_rowlen-1;
        Hgrid->d_nsw[0] = Hgrid->d_rowlen * Hgrid->d_numrows - 1;
        Hgrid->d_nse[0] = Hgrid->d_rowlen * (Hgrid->d_numrows-1) ;
        //bottom right hand corner
        Hgrid->d_nsw[Hgrid->d_rowlen-1] = Hgrid->d_rowlen*Hgrid->d_numrows - 1;
        Hgrid->d_nse[Hgrid->d_rowlen-1] = Hgrid->d_rowlen*Hgrid->d_numrows - 2;
        Hgrid->d_ne[Hgrid->d_rowlen-1] = 0;
        //rest of bottom row
        offset = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        std::cout << "bottom row offset " << offset << std::endl;
        for (unsigned int j=1; j<Hgrid->d_rowlen-1; j++) {
            Hgrid->d_nsw[j] = offset + j-1;
            Hgrid->d_nse[j] = offset + j;
        }
        //top row now
        //top left can resuse offset
        Hgrid->d_nsw[offset] = offset - 1;
        Hgrid->d_nw[offset] = offset + Hgrid->d_rowlen - 1;
        Hgrid->d_nnw[offset] = 0;
        Hgrid->d_nne[offset] = 1;
        //top right now
        offset = Hgrid->d_rowlen*Hgrid->d_numrows -1;
        Hgrid->d_nnw[offset] = Hgrid->d_rowlen-1;
        Hgrid->d_nne[offset] = 0;
        Hgrid->d_ne[offset] = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        //now the interior of the top row
        offset = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        std::cout << "top row offset " << offset << " d_nnw " << Hgrid->d_nnw[offset] << std::endl;
        std::cout << "bottom row 1 " << " d_nnw " << Hgrid->d_nnw[1] << std::endl;
        for (unsigned int j=1; j<Hgrid->d_rowlen-1; j++) {
            di = offset + j;
            Hgrid->d_nnw[di] = j;
            Hgrid->d_nne[di] = j+1;
        }
        std::cout << "bottom row 1 " << " d_nsw " << Hgrid->d_nsw[1] << std::endl;
        //now the b.c.s on the interior applied at start and end of each row
        for (unsigned int i=1; i<Hgrid->d_numrows-1; i+=2) {
            offset = i*Hgrid->d_rowlen;
            //left boundary
            Hgrid->d_nw[offset] = offset + Hgrid->d_rowlen - 1;
            //right boundary
            Hgrid->d_ne[offset+Hgrid->d_rowlen-1] = offset;
            Hgrid->d_nse[offset+Hgrid->d_rowlen-1] = offset - Hgrid->d_rowlen;
            Hgrid->d_nne[offset+Hgrid->d_rowlen-1] = offset + Hgrid->d_rowlen;
        }
        for (unsigned int i=2; i<Hgrid->d_numrows-2; i+=2) {
            offset = i*Hgrid->d_rowlen;
            //left boundary
            Hgrid->d_nnw[offset] = offset + 2*Hgrid->d_rowlen - 1;
            Hgrid->d_nw[offset] = offset + Hgrid->d_rowlen - 1;
            Hgrid->d_nsw[offset] = offset - Hgrid->d_rowlen - 1;
            //right boundary
            Hgrid->d_ne[offset+Hgrid->d_rowlen-1] = offset;
        }
    }

    //set boundary conditions with the bottom row even
    void setPeriodicEven() {

        int di = 0;
        int offset = Hgrid->d_rowlen*(Hgrid->d_numrows-1);
        FLT rowlen = Hgrid->d_rowlen;
        //along bottom row first
        //bottom left corner
        Hgrid->d_nw[di] = rowlen-1;
        Hgrid->d_nsw[di] = offset+rowlen-1;
        Hgrid->d_nse[di] = offset;
        //bottom right hand corner
        Hgrid->d_nsw[rowlen-1] = offset+rowlen-2;
        Hgrid->d_nse[rowlen-1] = offset+rowlen-1 ;
        Hgrid->d_ne[rowlen-1] = di;
        Hgrid->d_nne[rowlen-1] = rowlen;
        //rest of bottom row
        std::cout << "bottom row offset " << offset << std::endl;
        for (unsigned int j=1; j<Hgrid->d_rowlen-1; j++) {
            Hgrid->d_nsw[j] = offset + j-1;
            Hgrid->d_nse[j] = offset + j;
        }
        //top row now
        //top left can resuse offset
        Hgrid->d_nw[offset] = offset + rowlen -1;
        Hgrid->d_nnw[offset] = 0;
        Hgrid->d_nne[offset] = 1;
        //top right now
        Hgrid->d_nnw[offset+rowlen-1] = Hgrid->d_rowlen-1;
        Hgrid->d_nne[offset+rowlen-1] = 0;
        Hgrid->d_ne[offset+rowlen-1] = offset;
        Hgrid->d_nse[offset+rowlen-1] = offset-rowlen;
        //now the interior of the top row
        for (unsigned int j=1; j<Hgrid->d_rowlen-1; j++) {
            di = offset + j;
            Hgrid->d_nnw[di] = j;
            Hgrid->d_nne[di] = j+1;
        }
        std::cout << "top row offset " << offset << " d_nnw " << Hgrid->d_nnw[offset] << std::endl;
        std::cout << "bottom row 1 " << " d_nnw " << Hgrid->d_nnw[1] << std::endl;
        std::cout << "bottom row 1 " << " d_nsw " << Hgrid->d_nsw[1] << std::endl;
        //now the b.c.s on the interior applied at start and end of each row
        for (unsigned int i=2; i<Hgrid->d_numrows; i+=2) {
            offset = i*rowlen;
            std::cout << " left side offset " << offset << std::endl;
            //left boundary
            //Hgrid->d_nnw[offset] = offset + 2*Hgrid->d_rowlen - 1;
            Hgrid->d_nnw[offset] = offset + 2*rowlen - 1;
            Hgrid->d_nw[offset] = offset + rowlen - 1;
            Hgrid->d_nsw[offset] = offset - rowlen + 1;
            //right boundary
            Hgrid->d_ne[offset+rowlen-1] = offset;
        }
        for (unsigned int i=1; i<Hgrid->d_numrows; i+=2) {
            offset = i*rowlen;
            std::cout << " right side offset " << offset << std::endl;
            //left boundary
            Hgrid->d_nw[offset] = offset + rowlen - 1;
            //right boundary
            Hgrid->d_nne[offset+Hgrid->d_rowlen-1] = offset+rowlen;
            Hgrid->d_ne[offset+Hgrid->d_rowlen-1] = offset;
            Hgrid->d_nse[offset+Hgrid->d_rowlen-1] = offset-rowlen;
        }
    }

  //function to time step periodic b.c.s
    void step(FLT dt, FLT epsilon, FLT g, FLT k0, vector<complex<FLT>> oldPsi)
    {
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
        //cout << "in step before Stage 1" << endl;
        this->compute_dpsidt (this->psi, dpsidt, epsilon, g, k0);
        for (int di=0; di< this->n; ++di) {
            K1[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K1[di] * 0.5f ;
        }

        /*
         * Stage 2
        */
        //cout << "in step before Stage 2" << endl;
        this->compute_dpsidt (Ntst, dpsidt, epsilon, g, k0);
        for (int di=0; di< this->n; ++di) {
            K2[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K2[di] * 0.5f;
        }

        /*
         * Stage 3
         */
        //cout << "in step before Stage 3" << endl;
        this->compute_dpsidt (Ntst, dpsidt,  epsilon, g, k0);
        for (int di=0; di < this->n; ++di) {
            K3[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K3[di];
        }

        /*
         * Stage 4
         */
        //cout << "in step before Stage 4" << endl;
        this->compute_dpsidt (Ntst, dpsidt, epsilon, g, k0);
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

  //function to time step periodic b.c.s
    void stepHeat(FLT dt, FLT epsilon, FLT g, vector<complex<FLT>> oldPsi)
    {
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
        this->compute_dpsidtH (this->psi, dpsidt, epsilon, g);
        for (int di=0; di< this->n; ++di) {
            K1[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K1[di] * 0.5f ;
        }

        /*
         * Stage 2
        */
        this->compute_dpsidtH (Ntst, dpsidt, epsilon, g);
        for (int di=0; di< this->n; ++di) {
            K2[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K2[di] * 0.5f;
        }

        /*
         * Stage 3
         */
        this->compute_dpsidtH (Ntst, dpsidt, epsilon, g);
        for (int di=0; di < this->n; ++di) {
            K3[di] = dpsidt[di] * dt;
            Ntst[di] = this->psi[di] + K3[di];
        }

        /*
         * Stage 4
         */
        this->compute_dpsidtH (Ntst, dpsidt, epsilon, g);
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

    void check_green_boundary() {
        unsigned int offset = 0;
        for (unsigned int i=1; i<Hgrid->d_numrows-1; i++) {
            offset = i * Hgrid->d_rowlen;
            if (!(Hgrid->d_nnw[offset] == offset + 2*Hgrid->d_rowlen -1 && Hgrid->d_nw[offset] == offset + Hgrid->d_rowlen-1))
                std::cerr << "left boundary mismatch " << std::endl;;
            if (!(Hgrid->d_ne[offset+Hgrid->d_rowlen-1] == offset && Hgrid->d_nw[offset+Hgrid->d_rowlen-1] == offset - Hgrid->d_rowlen))
                std::cerr << "right boundary mismatch " << std::endl;
        }
    }


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


    vector<bool> complexZero(vector<FLT> invector, FLT min) {
        vector<bool> result;
        if (static_cast<int> (invector.size()) != this->n) {
            cerr << " compleZero invector " << invector.size() << " not equal to this->n " << this->n << endl;
            std::exit(0);
        }
        result.resize(invector.size());
        int count = 0;
        for(int di=0; di<this->n; di++) {
            if (invector[di] > min) {
                //std::cerr << "val " << invector[di] << " min " << min << std::endl;
                continue;
            }
            if (invector[di] > invector[this->Hgrid->d_ne[di]]) result[di] = false;
            else if (invector[di] > invector[this->Hgrid->d_nne[di]]) result[di] = false;
            else if (invector[di] > invector[this->Hgrid->d_nnw[di]]) result[di] = false;
            else if (invector[di] > invector[this->Hgrid->d_nw[di]]) result[di] = false;
            else if (invector[di] > invector[this->Hgrid->d_nsw[di]]) result[di] = false;
            else if (invector[di] > invector[this->Hgrid->d_nse[di]]) result[di] = false;
            else result[di] = true;
            if (result[di] == true) count++;
        }
        cout << " number of pinwheels " << count << " out of " << this->n << " hexes " << endl;
        return result;
    }


    vector<bool> contour(vector<FLT> invector) {
        vector<bool> result;
        if (static_cast<int> (invector.size()) != this->n) {
            cerr << " compleZero invector " << invector.size() << " not equal to this->n " << this->n << endl;
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
