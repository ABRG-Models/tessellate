#include <morph/tools.h>
#include <morph/Vector.h>
#include <morph/vVector.h>
//#include "opencv2/opencv.hpp"
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
#ifndef GENERAL
#include "topo/general.h"
#define GENERAL
#endif
#ifndef TOPO
#include "topo/topo.h"
#define TOPO
#endif
#include <morph/ShapeAnalysis.h>
// #include <boost/math/special_functions/bessel.hpp>

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::Tools;
using namespace std;
#ifndef PI
#define PI 3.14159265f
#endif

// to find the maximum value of a vector
class Analysis {
//global class variables
    public:
    struct extremum {
        int radialIndex;
        FLT radialValue;
    };
    vector<extremum> turnVal; //radial turning points
    vector<FLT> binVals;
    vector<FLT> histogram;
    vector<FLT> xs;
    vector<FLT> ys;
    int nBins, gaussBlur;
    CartHexSampler<FLT> C;

    FLT ROIwid, ROIpinwheelCount, patternFrequency, columnSpacing;

    //default constructor
    Analysis () {};

    //constructor for Fourier analysis
    Analysis (int nBins, int gaussBlur, FLT ROIwid) {
        this->nBins = nBins;
        this->gaussBlur = gaussBlur;
        this->ROIwid = ROIwid;
    };

    FLT maxVal( vector<FLT> invector) {
        FLT result = -1.0e7;
        for (unsigned int i = 0; i < invector.size(); i++) {
            if (invector[i] > result)
                result = invector[i];
        }
        return result;
    }


    // to find the minimum value of a vector
    FLT minVal( vector<FLT> invector) {
        FLT result = 1.0e7;
        for (unsigned int i = 0; i < invector.size(); i++) {
            if (invector[i] < result)
                result = invector[i];
            }
            return result;
    }


    // transform vector so its mean is zero
    vector<FLT> meanzero_vector(vector<FLT> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
        vector <FLT> result;
        int size = invector.size();
        //meanzero << "size " << size << endl;
        FLT sum = 0;
        FLT absSum = 0;
        for (int i=0; i <size; i++) {
            sum += invector[i];
            absSum += fabs(invector[i]);
        }
        sum = sum/(1.0*size);
        absSum = absSum / (1.0 * size);
        //meanzero << " mean  " << sum << endl;
        //meanzero << " absolute mean  " << absSum << endl;
        for (int i=0; i < size; i++) {
            result.push_back(invector[i] - sum);
            //meanzero << " i " << result[i] << endl;
        }
        return result;
    }

    // return the mean of a vector
        FLT mean_vector(vector<FLT> invector) {
        FLT result;
        int size = invector.size();
        //meanzero << "size " << size << endl;
        FLT sum = 0;
        for (int i=0; i <size; i++) {
            sum += invector[i];
        }
        sum = sum/(1.0*size);
        return result;
    }

    // return the mean of the absolute values of a  vector
        FLT absmean_vector(vector<FLT> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
          FLT result = 0;
	      FLT sum = 0;
          int size = invector.size();
        //meanzero << "size " << size << endl;
          for (int i=0; i <size; i++) {
            sum += fabs(invector[i]);
        }
        sum = sum/(1.0*size);
        return result;
    }
    //function to smooth a vector by moving average
    vector <FLT> smooth_vector(vector<FLT> invector, int window) {
        vector<FLT> outvector;
        int size = invector.size();
        outvector.resize(size);
        for (int i=1; i<size+1; i++) {
            outvector[i%size] = (invector[(i -1)%size] + invector[i%size] + invector[(i + 1)%size])/3.0;
        }
        return outvector;
    }

	vector <FLT> normalise (vector <FLT> invector) {
	  vector <FLT> result;
	  unsigned int size = invector.size();
	  result.resize(size, 0);
	  FLT maxV = -1e7;
	  FLT minV = 1e7;
	  for (unsigned int i=0;i<size;i++) {
	    if (invector[i] > maxV) {maxV = invector[i];}
		if (invector[i] < minV) {minV = invector[i];}
		}
	  FLT scaleV = 1./(maxV - minV);
	  for (unsigned int i=0;i<size;i++) {
	    result[i] = fmin(fmax((invector[i] - minV)*scaleV,0.),1.);
	  }
	  return result;
    }

/*
        vector <FLT> normalise (vector <FLT> invector) {
            vector <FLT> result;
            unsigned int size = invector.size();
            result.resize(size, 0);
            FLT maxV = -1e7;
            FLT minV = 1e7;
            for (unsigned int i=0;i<size;i++) {
                if (invector[i] > maxV) {maxV = invector[i];}
                    if (invector[i] < minV) {minV = invector[i];}
                }
            FLT scaleV = 1./(maxV - minV);
            for (unsigned int i=0;i<size;i++) {
                result[i] = fmin(fmax((invector[i] - minV)*scaleV,0.),1.);
            }
            return result;
        }
*/

    vector<FLT> getArg (vector<complex<FLT>> invVec) {
        vector<FLT> result;
        vector<complex<FLT>>::iterator p;
        for (p = invVec.begin(); p<invVec.end(); p++) {
            result.push_back(arg(*p));
        }
        return result;
    }


    vector<FLT> getArgPrincipal (vector<complex<FLT>> invVec) {
        vector<FLT> result;
        vector<complex<FLT>>::iterator p;
        for (p = invVec.begin(); p<invVec.end(); p++) {
            FLT phase;
            phase = arg(*p);
            if (phase >= 0.0f) {
                result.push_back(phase);
            }
            else {
                phase += 2.0f*PI;
                result.push_back(phase);
            }

        }
        return result;
    }


    vector<FLT> getArgPrincipalPi (vector<complex<FLT>> invVec) {
        vector<FLT> result;
        vector<complex<FLT>>::iterator p;
        for (p = invVec.begin(); p<invVec.end(); p++) {
            FLT phase;
            phase = arg(*p);
            if (phase >= 0.0f) {
                result.push_back(phase);
            }
            else {
                phase = phase + PI;
                result.push_back(phase);
            }

        }
        return result;
    }

    vector<FLT> getAbs (vector<complex<FLT>> invVec) {
        vector<FLT> result;
        vector<complex<FLT>>::iterator p;
        for (p = invVec.begin(); p<invVec.end(); p++) {
            result.push_back(abs(*p));
        }
        return result;
    }


    vector<FLT> getReal (vector<complex<FLT>> invVec) {
        vector<FLT> result;
        vector<complex<FLT>>::iterator p;
        for (p = invVec.begin(); p<invVec.end(); p++) {
            result.push_back(real(*p));
        }
        return result;
    }

    vector<FLT> getImag (vector<complex<FLT>> invVec) {
        vector<FLT> result;
        vector<complex<FLT>>::iterator p;
        for (p = invVec.begin(); p<invVec.end(); p++) {
            result.push_back(imag(*p));
        }
        return result;
    }

    vector<complex<FLT>> complexify(vector<FLT> r, vector<FLT> phase) {
        vector<complex<FLT>> result;
        result.resize(0);
        unsigned int rsize = r.size();
        //unsigned int psize = phase.size();
        result.resize(rsize);
        for (unsigned int i=0; i<rsize; i++) {
            std::complex temp = std::polar(r[i], phase[i]);
            result[i] = temp;
        }
        std::cout << "in complexify result size " << result.size() << " rsize " << rsize << std::endl;
        return result;
    }


  //function find_max to find turning points both values and indices.
    int find_max(vector<FLT> ray, int window) {
    int size = ray.size();
    // ofstream dfile ("turn.txt",ios::app);
    //dfile <<"size = " << size <<endl;
    vector<FLT> smoothRay;
    smoothRay = this->smooth_vector(ray, window);
    turnVal.resize(1000);
    //cout <<" "<<iend<<iend<<flush;
    FLT old_slope = 0;
    FLT new_slope = 0;
    int count = 0;
    old_slope = smoothRay[1] - smoothRay[0];
    for (int i =2; i<=size+1;i++){

      new_slope = smoothRay[i%size]-smoothRay[(i-1)%size];
      //dfile << " " << i%size << " " << old_slope << " "<<new_slope <<endl;
      if (new_slope*old_slope < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = smoothRay[i]; //should really interpolate
      //  dfile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
        count++;
      }
      old_slope = new_slope;
    }
    return count;
  }
// find the zeros in a ray angular
    int find_zeroDAngle(vector<int> ray) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
    int old_val = ray[0];
    int new_val;
    int count = 0;
    for (int i = 1 ; i<size+1;i++){
        if (ray[i%size] != 0) {
           new_val = ray[i%size];
           if (old_val*new_val == -1) {
                   count++;
           }
        old_val = new_val;
        }
      //zerofile << " radius " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }

// find the zeros in a ray angular
    vector<int> find_zeroIndices(vector<FLT> ray) {
    vector<int> result;
    result.resize(0);
    result.resize(0);
    int size = ray.size();
    ofstream zerofile ("zero.txt",ios::app);
    zerofile << "size = " << size << endl;
    FLT oldVal = ray[0];
    zerofile << " first oldVal = " << oldVal << " ray[0] " << ray[0] << endl;
    int count = 0;
    for (int i = 1 ; i<size+1; i++){
      //  if (ray[i%size] != 0.0) {
           FLT newVal = ray[i%size];
           zerofile << " radius " << i%size << " " << oldVal << " "<< newVal << " ray[i] " << ray[i%size] << endl;
           FLT norm = fabs(oldVal*newVal);
           if (oldVal*newVal/norm < 0.0) {
                 result.push_back(i);
                 count++;
           }
        oldVal = newVal;
        //}
      zerofile << " radius " << i%size << " " << oldVal << " "<< newVal << " ray[i] " << ray[i%size] << endl;
      }
      zerofile << "number of zeros " << count << endl;
      return result;
  }


 // find the zeros in a ray radial
    int find_zeroDRadius(vector<int> ray) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
    int old_val = ray[0];
    int new_val;
    int count = 0;
    for (int i = 1 ; i<size;i++){
        if (ray[i%size] != 0) {
           new_val = ray[i%size];
           if (old_val*new_val == -1) {
                   count++;
           }
        old_val = new_val;
        }
      //zerofile << " angle " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }


// find the zeros in a ray angular
    int find_zeroRadius(vector<int> ray) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
    int old_val = 0;
    int new_val = 0;
    int count = 0;
    for (int i =0; i<size;i++){
        if (ray[i%size] != 0) {
           new_val = ray[i%size];
           if (old_val*new_val == -1) {
                   count++;
           }
	   old_val = new_val;
        }
      //zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      }
    return count;
  }


    // find the zeros in a ray angular
    int find_zeroAngle(vector<FLT> ray, int window) {
    int size = ray.size();
    // ofstream zerofile ("zero.txt",ios::app);
    vector<FLT> smoothRay;
    //smoothRay = this->smooth_vector(ray, window);
    smoothRay = ray;
    //zerofile <<"size = " << size <<endl;
    turnVal.resize(1000);
    FLT old_val = 0;
    FLT new_val = 0;
    int count = 0;
    old_val = smoothRay[0];
    for (int i =1; i<size+1;i++){
      new_val = smoothRay[i%size];
      //zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      if (new_val*old_val < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue =smoothRay[i]; //should really interpolate
        //zerofile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
        count++;
      }
      old_val = new_val;
    }
    return count;
  }

      // find the zeros in a radial ray
    int find_zeroRadius(vector<FLT> ray, int window) {
    int size = ray.size();
    //ofstream zerofile ("zero.txt",ios::app);
     vector<FLT> smoothRay;
    // smoothRay = this->smooth_vector(ray, window);
     smoothRay = ray;
    //zerofile <<"size = " << size <<endl;

    turnVal.resize(1000);
    FLT old_val = 0;
    FLT new_val = 0;
    int count = 0;
    old_val = smoothRay[0];
    for (int i =1; i<size;i++){
      new_val = smoothRay[i%size];
      //zerofile << " " << i%size << " " << old_val << " "<<new_val <<endl;
      if (new_val*old_val < 0.) {
	turnVal[count].radialIndex = i;
	turnVal[count].radialValue = smoothRay[i]; //should really interpolate
        //zerofile << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
        count++;
      }
      old_val = new_val;
    }
    return count;
  }

    //scales a vector by a scalar
    vector<FLT> scaleVect(vector<FLT> invect, FLT scale) {
        vector<FLT> result;
        int size = invect.size();
        result.resize(size);
        for (int i=0; i<size; i++) {
            result[i] = invect[i]*scale;
        }
        return result;
    }

    //counts true vals in bool array
    int countBool(std::vector<bool> inVect) {
        int count = 0;
        for (unsigned int i=0; i<inVect.size(); i++) {
            if (inVect[i]) count++;
        }
        return count;
    }

    vector <FLT> lengthenVector (vector<FLT> svector, int lSize) {
        vector <FLT> result;
        int sSize = 0;
        FLT sStep, lStep = 0;
        result.resize(0);
        sSize = svector.size();
        sStep = 1.0 / (1.0 * (sSize-1));
        lStep = 1.0 / (1.0 * (lSize-1));
        FLT start = 0;
        FLT finish = 0;
        FLT value = 0;
        FLT delta = 0.0000001;
        int marker = 0;
        for (int i=0; i<sSize-1; i++) { // walk along the short vector
            start = i*sStep;
            finish = (i+1)*sStep + delta;
            while ((marker  < lSize) && (marker*lStep < finish)) { //walk along the long vector
                    value = (svector[i+1]*(marker*lStep - start) + svector[i]*(finish - marker*lStep))/sStep;
                    result.push_back(value);
                    marker++;
            }
        }
        if (marker != lSize){
            std::cout <<  " lSize " << lSize << " sSize " << sSize << " count " << marker <<  " not filled" << std::endl;
            result.resize(0);
            return result;
        }
        else {
            //std::cout << " l size " << lSize << " sSize " << sSize << "resSize " << result.size() << endl;
            // result.push_back(svector[sSize - 1]);
            return result;
        }
    } //end of function lengthenVector


    //Convert a vect to cv::Mat for a rectangular grid
    cv::Mat vect2mat(std::vector<FLT> A, int nrows, int rowl) {

        //C.nx = nx;
        //C.ny = ny;
        int n = nrows*rowl;
        FLT maxV = -1e9;
        FLT minV = +1e9;
        for(int i=0; i<n; i++){
            if(maxV<fabs(A[i])){
                maxV = fabs(A[i]);
            } // THIS WAY ENSURE THAT ZERO DIFF ALWAYS MAPS TO VAL OF 0.5
            if(minV>A[i]){ minV = A[i]; }
        }
        FLT scale = 1./(2*maxV);

        cv::Mat I = cv::Mat::zeros(nrows,nrows,CV_32F);
        vector<FLT> row;
        vector<FLT> longrow;
        row.resize(rowl);
        longrow.resize(nrows);
        int k=0;
        for(int i=0;i<nrows;i++){
            for(int j=0; j<rowl; j++){
                row[j] = (A[k]+maxV)*scale;
                k++;
            }
            longrow = lengthenVector(row, nrows);
            for (int j=0; j<nrows; j++) {
                I.at<FLT>(i,j) = longrow[j];
            }
        }
        std::cout << "in vect2mat k " << k << std::endl;
        return I;
    }


    //Convert a vect to cv::Mat for a parallelogram grid
    //this will convert the parallelogram to a rectangle
    cv::Mat vect2matPar(std::vector<FLT> A, int nrows, int rowl) {

        int n = nrows*rowl;
        FLT maxV = -1e9;
        FLT minV = +1e9;
        for(int i=0; i<n; i++){
            if(maxV<fabs(A[i])){
                maxV = fabs(A[i]);
            }
            if(minV>A[i]){ minV = A[i]; }
        }
        FLT scale = 1./(1.0*maxV);
        vector<FLT> col;
        FLT pgram[nrows][nrows];
        cv::Mat I = cv::Mat::zeros(nrows+2,nrows+2,CV_32F);
        int k=0;
        std::cout << "in vect2matPar nrows " << nrows << " rowl " << " scale " << scale << std::endl;
        for(int j=0;j<nrows;j++){
            for(int i=0; i<nrows; i++){
                int offset = i*rowl;
                pgram[i][j] = (A[offset+j])*scale;
                k++;
            }

        }

        int offset = (nrows-1)/2 + 1;
        for(int i=0;i<nrows;i+=2){
            offset--;
            for(int j=0; j<nrows; j++) {
                I.at<FLT>(i,j) = pgram[i][(j+offset)%rowl];
            }
            int ir = i+1;
            for(int j=0; j<nrows; j++){
                I.at<FLT>(ir,j) = pgram[ir][(j+offset)%rowl];
            }
            //std::cout <<"in vect2matPar i " << i <<  " offset " << offset << " k= " << k <<std::endl;
        }
        return I;
    }

    //Convert a vect to cv::Mat for a parallelogram grid
    //cutting out a square region on the base.
    cv::Mat vect2matCut(std::vector<FLT> A, int nrows, int rowl) {

        int n = nrows*rowl;
        FLT maxV = -1e9;
        FLT minV = +1e9;
        for(int i=0; i<n; i++){
            if(maxV<fabs(A[i])){
                maxV = fabs(A[i]);
            }
            if(minV>A[i]){ minV = A[i]; }
        }
        FLT scale = 1./(1.0*maxV);
        vector<FLT> col;
        col.resize(nrows);
        FLT pgram[nrows][rowl];
        int k=0;
        //std::cout << "in vect2matPar nrows " << nrows << " rowl " << rowl << " longnrows " << longnrows << std::endl;
        for(int j=0;j<rowl;j++){
            for(int i=0; i<nrows; i++){
                int offset = i*rowl;
                pgram[i][j] = (A[offset+j])*scale;
                k++;
            }

        }

        int offset = (nrows-1)/2 + 1;
        int trail = 0;
        //cv::Mat I = cv::Mat::zeros(longnrows,rowl,CV_32F);
        cv::Mat I = cv::Mat::zeros(rowl-offset+2,rowl-offset+2,CV_32F);
        for(int i=0;i<rowl-(nrows-1)/2;i+=2){
            offset--;
            trail++;
            for(int j=offset; j<rowl-trail; j++) {
                I.at<FLT>(i,j-offset) = pgram[i][j];
            }
            int ir = i+1;
            for(int j=offset; j<rowl-trail; j++){
                I.at<FLT>(ir,j-offset) = pgram[ir][j];
            }
            //std::cout <<"in vect2matCut i " << i <<  " offset " << offset << "trail =  " << trail  <<std::endl;
        }
        return I;
    }

    cv::Mat sqmatrix2mat(vector<vector<FLT>> matrix, int iRoI) {
        std::cout << "in sqmatrix2mat size " << iRoI << " matrix size " << matrix.size() << std::endl;
        cv::Mat I = cv::Mat::zeros(iRoI, iRoI, CV_32F);
        for (int i=0; i<iRoI; i++) {
            for (int j=0; j<iRoI; j++) {
                I.at<FLT>(i,j) = matrix[i][j];
               // std::cout << "ij" << i << j << " matrix val " << matrix[i][j] << " I val " << I.at<FLT>(i,j) << std::endl;
            }
        }
        return I;
    }

    //Convert a vect to a vector<vector<FLT>> for a parallelogram grid
    vector<vector<FLT>> vect2img(std::vector<FLT> A, int nrows, int rowl) {
        vector<vector<FLT>> result;
        int n = nrows*rowl;
        FLT maxV = -1e9;
        FLT minV = +1e9;
        for(int i=0; i<n; i++){
            if(maxV<fabs(A[i])){
                maxV = fabs(A[i]);
            }
            if(minV>A[i]){ minV = A[i]; }
        }
        FLT scale = 1./(1.0*maxV);
        int longnrows = floor(1.29*nrows);
        //C.nx = rowl;
        //C.ny = longnrows;
        vector<FLT> col;
        vector<FLT> longcol;
        col.resize(nrows);
        longcol.resize(longnrows);
        result.resize(longnrows);
        FLT pgram[longnrows][rowl];
        int k=0;
        //std::cout << "in vect2img nrows " << nrows << " rowl " << rowl << " longnrows " << longnrows << std::endl;
        for(int j=0;j<rowl;j++){
            for(int i=0; i<nrows; i++){
                int offset = i*rowl;
                col[i] = (A[offset+j])*scale;
                k++;
            }
            longcol = lengthenVector(col,longnrows);
            for (int i=0; i<longnrows; i++) {
                pgram[i][j] = longcol[i];
            }
        }

        int offset = (longnrows-1)/2 + 1;
        for(int i=0;i<longnrows;i+=2){
            offset--;
            for(int j=0; j<rowl; j++) {
                result[i].push_back(pgram[i][(j+offset)%rowl]);
            }
            int ir = i+1;
            for(int j=0; j<rowl; j++){
                result[i].push_back(pgram[ir][(j+offset)%rowl]);
            }
            //std::cout <<"in vect2matPar i " << i <<  " offset " << offset << " k= " << k <<std::endl;
        }
        return result;
    }

    std::vector<FLT> getPatternFrequency(cv::Mat I, bool showfft) {
        // ANALYSIS STEP 5. ESTIMATE ISO-ORIENTATION COLUMN SPACING
        this->histogram.resize(this->nBins);
        this->binVals.resize(this->nBins,0.f);
        std::cout << "in getpatternFrequency nbins " << this->nBins << std::endl;

        std::vector<std::vector<FLT> > h1 = C.fft(I, this->nBins, this->gaussBlur, showfft);

        //std::cout << "in getpatternFrequency " << std::endl;
        binVals = h1[0];      // get histogram bin mid-values
        histogram = h1[1];
        // sample portion of histogram to fit
        int nsamp = this->nBins;
        int ioff = 1;
        this->xs.resize(nsamp,0.0);
        this->ys.resize(nsamp,0.0);
        arma::vec Xs(nsamp-ioff);
        arma::vec Ys(nsamp-ioff);
        for(int i=0; i<nsamp; i++){
            xs[i] = binVals[i];
            ys[i] = histogram[i];
            //std::cout << "xs " << binVals[i] << " ys " << histogram[i] << std::endl;
        }

        for (int i=0; i<nsamp-ioff; i++) {
            Xs[i] = binVals[i+ioff];
            Ys[i] = histogram[i+ioff];
        }

        int polyOrder = 10;
        // do polynomial fit
        arma::vec cf = arma::polyfit(Xs,Ys,polyOrder);

        // make a high-resolution model for the data
        arma::vec xfit(nsamp);
        for(int i=0;i<nsamp;i++){
            xfit[i] = binVals[nsamp-1]*i/(1.0*(nsamp-1));
        }
        arma::vec yfit = arma::polyval(cf,xfit);
/*
        // get frequency at which high-res model peaks
        FLT maxVal = -1e9;
        FLT maxX = 0;
//this is cludged in order to avoid the zero frequency peak.
        for(int i=3;i<nsamp;i++){
            if(ys[i]>maxVal){
                maxVal = ys[i];
                maxX = xs[i];
            }
        }
*/


        // get frequency at which high-res model peaks
        float maxVal = -1e9;
        float maxX = 0;
        for(int i=1;i<nsamp;i++){
            //std::cout << "yfit " << i << " is " << yfit[i] << " xfit " << xfit[i] << std::endl;
            if(yfit[i]>maxVal){
                maxVal = yfit[i];
                maxX = xfit[i];
            }
        }


        // return coeffs in standard vector
        std::vector<float> coeffs(cf.size());
        for(int i=0;i<cf.size();i++){
            coeffs[i] = 1.0f*cf[i];
        }


        this->patternFrequency = maxX; // units are cycles / ROI-width
        this->columnSpacing = this->ROIwid / patternFrequency;  // spacing between iso-orientation columns in units of cortex sheet, e.g., to plot scale bar on maps
        std::cout <<"in getPatternFrequency MaxVal " << maxVal << " maxX " << maxX << std::endl;
        return coeffs;
    }


  //function bessel_ray for computing bessel functions along a ray
  // takes a vector of FLTs representing the radius
  // returns the Bessel function values
  // vector<FLT> bessel_ray (int v, vector<FLT> ray) {
  //   vector <FLT> result(n,0.);
  //   result = boost::cyl_bessel_j(v, ray);
  //   return result;
  // }
}; //end of class Analysis
