#include <morph/tools.h>
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

// to find the maximum value of a vector
class Analysis {
//global class variables
public:
  struct extremum {
  int radialIndex;
  FLT radialValue;
  };
  vector<extremum> turnVal; //radial turning points

  //default constructor

  Analysis () {};

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
/*!
 * for mixing bits of three arguments used to generate a good random
 * seed using time(), getpid() and clock()
 */
unsigned int
mix (unsigned int a, unsigned int b, unsigned int c)
{
    a=a-b;  a=a-c;  a=a^(c >> 13);
    b=b-c;  b=b-a;  b=b^(a << 8);
    c=c-a;  c=c-b;  c=c^(b >> 13);
    a=a-b;  a=a-c;  a=a^(c >> 12);
    b=b-c;  b=b-a;  b=b^(a << 16);
    c=c-a;  c=c-b;  c=c^(b >> 5);
    a=a-b;  a=a-c;  a=a^(c >> 3);
    b=b-c;  b=b-a;  b=b^(a << 10);
    c=c-a;  c=c-b;  c=c^(b >> 15);
    return c;
}



  //function find_max to find turning points both values and indices.
    int find_extrema(vector<FLT> ray, int window=0) {
        int size = ray.size();
        vector<FLT> smoothRay;
        smoothRay = this->smooth_vector(ray, window);
        turnVal.resize(1000);
        FLT old_slope = 0;
        FLT new_slope = 0;
        int count = 0;
        old_slope = smoothRay[1] - smoothRay[0];
        for (int i =2; i<=size+1;i++){
            new_slope = smoothRay[i%size]-smoothRay[(i-1)%size];
            if (new_slope*old_slope < 0.) {
                turnVal[count].radialIndex = i;
                turnVal[count].radialValue = smoothRay[i]; //should really interpolate
                std::cout << "turn index " << turnVal[count].radialIndex << " turn value " << turnVal[count].radialValue << endl;
                count++;
            }
            old_slope = new_slope;
        }
        return count;
    }

  //function return indices of maxima in an array.
    vector<int> find_maxIndices(vector<FLT> ray, int window=0) {
    vector<int> result;
    result.resize(0);
    int size = ray.size();
    vector<FLT> smoothRay;
    smoothRay = this->smooth_vector(ray, window);
    FLT old_slope = 0.0;
    FLT new_slope = 0.0;
    old_slope = smoothRay[1] - smoothRay[0];
    for (int i=2; i<=size+2; i++){
        new_slope = smoothRay[i%size]-smoothRay[(i-1)%size];
        if (old_slope < 0.0 && new_slope > 0.0) {
            result.push_back(i);
            std::cout << "turn index " << i << " turn value " << (smoothRay[i%size] + smoothRay[(i+1)%size])/2.0 << " old slope " << old_slope << " new slope " << new_slope  << endl;
        }
      old_slope = new_slope;
    }
    return result;;
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

// find the indices of zeros in a ray angular
    vector<int> find_zeroIndices(vector<FLT> ray) {
    vector<int> result;
    result.resize(0);
    result.resize(0);
    int size = ray.size();
    ofstream zerofile ("./zero.txt",ios::app);
    zerofile << "size = " << size << endl;
    FLT oldVal = ray[0];
    zerofile << " first oldVal = " << oldVal << " ray[0] " << ray[0] << endl;
    int count = 0;
    for (int i = 1 ; i<size+1; i++){
      //  if (ray[i%size] != 0.0) {
           FLT newVal = ray[i%size];
           FLT norm = fabs(oldVal*newVal);
           if (oldVal*newVal/norm < 0.0 && oldVal <0.0) {
                zerofile << " radius " << i%size << " " << oldVal << " "<< newVal << " ray[i] " << ray[i%size] << endl;
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
      if (new_val*old_val < 0.0) {
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

    //find the Euclidean norm of a vector
    FLT vectNorm(std::vector<FLT> invector)
    {
        FLT result = 0;
        for (unsigned int i=0; i< invector.size(); i++) {
            result += invector[i]*invector[i];
        }
        return sqrt(result);
    }



    //find the normed difference between two vectors
    //TOFIX this assumes that the vectors are the same size, need error handling
    FLT  normedDiff(vector<FLT> preField, vector<FLT> currField)
    {
        FLT result;
        FLT diffField = 0;
        unsigned int fieldSize = preField.size();
        for (unsigned int i=0; i<fieldSize; i++) {
            diffField += ((currField[i] - preField[i]) * (currField[i] - preField[i]));
        }
        FLT preNorm = this->vectNorm(preField);
        FLT currNorm = this->vectNorm(currField);
        if (preNorm * currNorm > 0.0) {
            result = diffField / (currNorm * currNorm);
            return result;
        }
        else {
            result = 0;
            return result;
        }
    }//end of method normedDiff


  //function bessel_ray for computing bessel functions along a ray
  // takes a vector of FLTs representing the radius
  // returns the Bessel function values
  // vector<FLT> bessel_ray (int v, vector<FLT> ray) {
  //   vector <FLT> result(n,0.);
  //   result = boost::cyl_bessel_j(v, ray);
  //   return result;
  // }
}; //end of class Analysis
