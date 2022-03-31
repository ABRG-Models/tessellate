#include <morph/tools.h>
#include <morph/Vector.h>
#include <morph/vVector.h>
#include "opencv2/opencv.hpp"
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
    vector<FLT> binVals,histogram;
    int nBins, gaussBlur;
    CartHexSampler<FLT> CHM;
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

    vector<FLT> getAbs (vector<complex<FLT>> invVec) {
        vector<FLT> result;
        vector<complex<FLT>>::iterator p;
        for (p = invVec.begin(); p<invVec.end(); p++) {
            result.push_back(abs(*p));
        }
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

    //counts true vals in bool array
    int countBool(std::vector<bool> inVect) {
        int count = 0;
        for (unsigned int i=0; i<inVect.size(); i++) {
            if (inVect[i]) count++;
            return count;
        }
    }


    //Convert a vect to cv::Mat
    cv::Mat vect2mat(std::vector<FLT> A, int nx, int ny) {

        CHM.C.nx = nx;
        CHM.C.ny = ny;
        int n = nx*ny;
        FLT maxV = -1e9;
        FLT minV = +1e9;
        for(int i=0; i<n; i++){
            if(maxV<fabs(A[i])){
                maxV = fabs(A[i]);
            } // THIS WAY ENSURE THAT ZERO DIFF ALWAYS MAPS TO VAL OF 0.5
            if(minV>A[i]){ minV = A[i]; }
        }
        FLT scale = 1./(2*maxV);

        cv::Mat I = cv::Mat::zeros(nx,ny,CV_32F);
        int k=0;
        for(int i=0;i<nx;i++){
            for(int j=ny; j-->0;){ // invert y axis of image
                I.at<float>(j,i) = (A[k]+maxV)*scale;
                k++;
            }
        }
        return I;
    }

    std::vector<FLT> getPatternFrequency(cv::Mat I, bool showfft){

        // ANALYSIS STEP 5. ESTIMATE ISO-ORIENTATION COLUMN SPACING
        int sampleRange = 1;
        int polyOrder = 4;
        binVals.resize(nBins,0.);
        histogram.resize(nBins,0.);

        // Get frequency histogram from image
        //cv::Mat I1 = CHM.getDifferenceImage(orResponseSampled[0],orResponseSampled[2]);
        std::vector<std::vector<FLT> > h1 = CHM.fft(I, nBins, gaussBlur, showfft);

        // sample portion of histogram to fit
        int nsamp = nBins*sampleRange;
        arma::vec xs(nsamp);
        arma::vec ys(nsamp);
        for(int i=0;i<nsamp;i++){
            xs[i] = binVals[i];
            ys[i] = histogram[i];
        }

        // do polynomial fit
        arma::vec cf = arma::polyfit(xs,ys,polyOrder);

        // make a high-resolution model for the data
        int fitres = 1000;
        arma::vec xfit(fitres);
        for(int i=0;i<fitres;i++){
            xfit[i] = binVals[nsamp-1]*(FLT)i/((FLT)fitres-1);
        }
        arma::vec yfit = arma::polyval(cf,xfit);

        // get frequency at which high-res model peaks
        FLT maxVal = -1e9;
        FLT maxX = 0;
        for(int i=0;i<fitres;i++){
            if(yfit[i]>maxVal){
                maxVal = yfit[i];
                maxX = xfit[i];
            }
        }

        this->patternFrequency = maxX; // units are cycles / ROI-width
        this->columnSpacing = ROIwid / patternFrequency;  // spacing between iso-orientation columns in units of cortex sheet, e.g., to plot scale bar on maps

        // return coeffs in standard vector
        std::vector<FLT> coeffs(cf.size());
        for(int i=0;i<cf.size();i++){
            coeffs[i] = (FLT)cf[i];
        }

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
