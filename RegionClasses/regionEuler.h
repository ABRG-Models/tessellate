/* 
 * Dregion class
 * Author: John Brooke
 *
 * Date 2019/10
 *
 * Creates a hexgrid, divides it into Dirichlet domains
 * then solves KS equations and can morph to curved 
 * boundaries
 *
 */
#include <morph/tools.h>
#include <morph/HexGrid.h>
#include <morph/ReadCurves.h>
#include <morph/HdfData.h>
#include <morph/BezCurve.h>
#include <morph/BezCurvePath.h>
#include <morph/BezCoord.h>
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
#include "hexGeometry.h"
// #include <boost/math/special_functions/bessel.hpp>
#define PI 3.1415926535897932
//#define NUMPOINTS 5 //just the A-E rows.
#define NUMPOINTS 79 //just the A-E rows.

using std::vector;
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
using morph::HexGrid;
using morph::HdfData;
using morph::ReadCurves;
using morph::Tools;
using morph::BezCurve;
using morph::BezCurvePath;
using morph::BezCoord;
using namespace std;

#define NE(hi) (this->Hgrid->d_ne[hi])
#define HAS_NE(hi) (this->Hgrid->d_ne[hi] == -1 ? false : true)

#define NW(hi) (this->Hgrid->d_nw[hi])
#define HAS_NW(hi) (this->Hgrid->d_nw[hi] == -1 ? false : true)

#define NNE(hi) (this->Hgrid->d_nne[hi])
#define HAS_NNE(hi) (this->Hgrid->d_nne[hi] == -1 ? false : true)

#define NNW(hi) (this->Hgrid->d_nnw[hi])
#define HAS_NNW(hi) (this->Hgrid->d_nnw[hi] == -1 ? false : true)

#define NSE(hi) (this->Hgrid->d_nse[hi])
#define HAS_NSE(hi) (this->Hgrid->d_nse[hi] == -1 ? false : true)

#define NSW(hi) (this->Hgrid->d_nsw[hi])
#define HAS_NSW(hi) (this->Hgrid->d_nsw[hi] == -1 ? false : true)
  
class DRegion
{
public:
    int scale;
    int n;
    double ds;
	//double overds;
  struct point {
    double xval;
    double yval;
  };

  // list of objects visible to member functions
  //vector<double> rad;
  vector<vector<int> > N; // hex neighbourhood 
  vector<vector<int> > region; //for each hex sorted list of regions
  vector<list<Hex>> regionIndex; //for each region list of hexes it contains
  vector<vector<int>> hexRegionList; //for each hex neighbour regions
  vector<vector<int>> regionList; //for each region, regions that are its neighbours
  vector<vector<int>> regionVertex; //for each region, vertices that bound it
  vector<vector<int>> sortedBoundary; // indices of boundary hexes sorted by angle
  vector<int> edgeIndex; // vector of the keys for the integers representing the edge pairs
  vector<list<Hex>> regionBound; //for each region, index of boundary vertices
  std::map<int,vector<int>> edges; //map of (i,j) edges, uses pair<->int converters
  vector<vector<double> > regionDist; //from each hexdistances to each seed point
  vector<vector<pair<double,double>>> vCoords;
  vector<vector<pair<double,double>>> mCoords;
  //vector<double> psi;
  vector<int> Creg; //for each count of regions it touches
  //vector<int> Cnbr; //for each hex count of neighbour hexes
  vector<double> NN, CC; //hold the field values for each hex
  pair <float, float> centres[NUMPOINTS]; //seed points for regions
  vector<std::pair<double,double>> diff; //difference between seed point and CoG of region
  morph::HexGrid* Hgrid;
  vector<morph::BezCurvePath<float>> curvedBoundary;
  hexGeometry* hGeo;
  vector<vector<hexGeometry::lineSegment>> radialSegments; //radial segements from original vertices

  int base = 1000;
  //class constructor
  DRegion (int scale, string logpath) {
  ofstream afile (logpath + "/debug.out" );
  ofstream jfile (logpath + "/correlateVector.out");
  cout << "before creating BezCurve" <<endl;;
  srand(time(NULL));
    this->scale = scale;
    double s = pow(2.0, scale-1);
	ds = 1.0/s;
  n = 0;
  Hgrid = new HexGrid(this->ds, 2.0, 0.0, morph::HexDomainShape::Boundary);
  n = Hgrid->num();
  cout << "after creating HexGrid"<<endl;
  //double maxX = Hgrid->getXmax(0.0);
  //double minX = Hgrid->getXmin(0.0);
  cout << "before filling H " << Hgrid->num() << endl;
  hGeo = new hexGeometry();
// now read in the boundary either as a header or as a morph read
#include "bezRectangle.h"
//  morph::ReadCurves r("./barrelAE.svg");
//  Hgrid->setBoundary (r.getCorticalPath());
// this was the original call, I am trying out setBoundaryDregion for debugging 
    Hgrid->setBoundaryDRegion(bound);
//  Hgrid->setBoundary (bound);
  cout << "after setting boundary on  H " << Hgrid->num() << endl;
  n = Hgrid->num();
  cout << "after  filling H " << " n = " << n <<endl;
// now set the centres either read in or randomly generated
  #include "centres.h"
  /*
  for (int i=0;i<NUMPOINTS;i++) {
      double choice = morph::Tools::randDouble();
	  if ((0 < choice) && (choice <= 0.25)) {
          centres[i].first = (morph::Tools::randDouble()) * maxX;
          centres[i].second = (morph::Tools::randDouble()) * maxY  ;
	  }	
      else if ((0.25 < choice) && (choice <= 0.5))
	  {
          centres[i].first = (morph::Tools::randDouble()) * maxX;
          centres[i].second = (morph::Tools::randDouble()) * minY  ;
	  }
	  else if ((0.5 < choice) && (choice  <= 0.75)) {
          centres[i].first = (morph::Tools::randDouble()) * minX;
          centres[i].second = (morph::Tools::randDouble()) * maxY  ;
	  }	
      else 
	  {
          centres[i].first = (morph::Tools::randDouble()) * minX;
          centres[i].second = (morph::Tools::randDouble()) * minY;
	  }
   } //end of setting of random values
*/
   cout << "after setting centres " << endl;
//these are the vectors of vectors for the regions
  regionDist.resize(n);
  region.resize(n);  
  N.resize(n);
  curvedBoundary.resize(NUMPOINTS);
  cout << "after  filleting fish " << " n = " << n <<endl;
  //this->Cnbr.resize(n,6); //count of neighbouring hexes
  this->Creg.resize(n,0); //count of neighbouring hexes
  cout << "before hexRegionList"<<endl;
  hexRegionList.resize(n); //neighbouring regions for a hex
  regionList.resize(NUMPOINTS); //neighbouring regions for a region
  regionVertex.resize(NUMPOINTS); //vertices for a region
  radialSegments.resize(NUMPOINTS);
  regionBound.resize(NUMPOINTS);
  sortedBoundary.resize(NUMPOINTS);
  vCoords.resize(NUMPOINTS);
  mCoords.resize(NUMPOINTS);
  cout << "before neighbour array" << endl;  
// making a neighbour array for convenience
  for (int idx = 0; idx < n; idx++) {
    N[idx].resize(6);
    N[idx][0] = Hgrid->d_ne[idx];
    N[idx][1] = Hgrid->d_nne[idx];
    N[idx][2] = Hgrid->d_nnw[idx];
    N[idx][3] = Hgrid->d_nw[idx];
    N[idx][4] = Hgrid->d_nsw[idx];
    N[idx][5] = Hgrid->d_nse[idx];
  }
cout << "after neighbour array" << endl;  

	//arrays to hold the field values
    NN.resize(n);
    CC.resize(n);
	cout << "after alloc NN and CC" <<endl;

    vector <vector <double> > sortedDist;
    sortedDist.resize(n);
	// get a list of distances from each Dirichlet point for each hex.
    for (auto h : Hgrid->hexen) {
      for (int j=0;j<NUMPOINTS;j++) {
        double temp = h.distanceFrom(centres[j]);
    //cout << "after access h " << h.vi  <<endl;
        // afile << "   " <<hexcen.xval <<"  "<<hexcen.yval;
        regionDist[h.vi].push_back(temp);
        // afile << temp <<"  ";
      }
    }
    //cout << "after regionDist"  << endl;
    // to produce a list of the sorted distances of each hex from the seed points
    // Note the process of sorting the indexes in distance order is completely independent
    // from the production of Sorted list but we assume that the similar sorting used for
    // each means that they are compatible.

    for (int i=0;i<n;i++) {
     	vector <double> tempvector1;
        tempvector1 = regionDist[i];
     	std::stable_sort(tempvector1.begin(),tempvector1.end());
     	sortedDist[i] = tempvector1;
        vector <int> tempint = sort_indexes(regionDist[i]);
        region[i] = tempint;

      }
    cout << "after region[i] allocated " <<endl;

     //afile <<"list of all the regions"<<endl;
     for(int i=0;i<n;i++){

       if (sortedDist[i][0] == sortedDist[i][1]){
          afile<<"i ="<<i;
     	  afile << "  "<<region[i][0]<<" is equidistant from "<< region[i][1];
       }
     }

     cout << "before hex list loop" << endl;
     // this determines the type of hex
      for (auto h : Hgrid->hexen){
	    //Cnbr.at(h.vi) = 6;

    //      if(h.boundaryHex == true) {
		    //cout << " its a boundary hex i= " << h.vi << endl;
            if (!HAS_NE(h.vi)) {
				hexRegionList[h.vi].push_back(-1);
		        h.setBoundaryHex();
				cout << "no ne" << endl;
            } 
            else {
              hexRegionList[h.vi].push_back(region[N[h.vi][0]][0]);//nbr region 
            }
                  
            if (!HAS_NNE(h.vi)) {
				hexRegionList[h.vi].push_back(-1);
		        h.setBoundaryHex();
				cout << "no nne" << endl;
            } 
            else {
              hexRegionList[h.vi].push_back(region[N[h.vi][1]][0]); //nbr region 
            }
            if (!HAS_NNW(h.vi)) {
				hexRegionList[h.vi].push_back(-1);
				cout << "no nnw" << endl;
		        h.setBoundaryHex();
            } 
            else {
              hexRegionList[h.vi].push_back(region[N[h.vi][2]][0]); //nbr region
            }
            if (!HAS_NW(h.vi)) {
				hexRegionList[h.vi].push_back(-1);
		        h.setBoundaryHex();
				cout << "no nw" << endl;
            } 
            else {
              hexRegionList[h.vi].push_back(region[N[h.vi][3]][0]); //nbr region 
            }
            if (!HAS_NSW(h.vi)) {
				hexRegionList[h.vi].push_back(-1);
		        h.setBoundaryHex();
				cout << "no nsw" << endl;
            } 
            else {
              hexRegionList[h.vi].push_back(region[N[h.vi][4]][0]); //nbr region
            }
            if (!HAS_NSE(h.vi)) {
				hexRegionList[h.vi].push_back(-1);
				h.setBoundaryHex();
				cout << "no nse" << endl;
            } 
            else {
              hexRegionList[h.vi].push_back(region[N[h.vi][5]][0]); //nbr region 
            }
			/*
		  }
          else {
		       // cout << "its not a boundary hex" << endl;
               for(int j=0;j<6;j++) {
                  hexRegionList[h.vi].push_back(region[N[h.vi][j]][0]); //nbr region
               }
      //       } //end of if dealing with outer boundary
	      */
		  // processing of internal boundaries
          int centralRegion = region[h.vi][0];
          int oldRegion = centralRegion;
          int newRegion = 0;
 cout << "just before the internal boundary logic" << " i " << h.vi << endl;
          for(int j=0;j<6;j++) {
		     //cout << "j = " << j  << " h.vi " << h.vi <<  endl;
             newRegion =  hexRegionList[h.vi][j];
		     cout << " hexRegionList " << hexRegionList[h.vi][j] << endl;
             if (centralRegion != newRegion){ //its a boundary hex
			   cout << "centralRegion " << centralRegion << " newRegion " << newRegion << endl;
               h.setBoundaryHex();
			   N[h.vi][j] = h.vi;
               //cout << " Cnbr " << Cnbr[h.vi] << endl;
               if (oldRegion != newRegion){ //logic to test if vertex
                 Creg[h.vi]++;
				 cout << " Creg " << Creg[h.vi] << " oldRegion " << oldRegion << " newRegion " << newRegion << endl;
                 oldRegion = newRegion;
               }
             }
           }
		    cout << "end of creg loop" << endl;
      } //end of logic determining hex types

		this->regionIndex.resize(NUMPOINTS);
		for (int j=0;j<NUMPOINTS;j++){
        for (auto h : Hgrid->hexen) {
           //count++;
           //cout << "hex " << count << endl;
	       if (region[h.vi][0] == j) {
           cout << "hex " << h.vi << " region " << region[h.vi][0] << " j " << j<< endl;
      	   this->regionIndex[j].push_back(h);
        }
	  }
	}


      cout << "after creating internal boundaries" << endl;
/*
	  //now mark all boundary hexes, internal or external as onBoundary
	  for (auto h : Hgrid->hexen) {
	    if ((Cnbr[h.vi] < 6) || (Creg[h.vi] == 1)) {
		  h.setBoundaryHex();
		  }
	  }
*/



      diff.resize(NUMPOINTS);
	  int totalHex=0; 
      for (int j=0;j<NUMPOINTS;j++){
	      totalHex += printRegion(j); 
		  cout<<endl;
		  }
      for (int j=0;j<NUMPOINTS;j++){
          diff[j] = this->set_polars(j);
          afile << "diff seed-centre" << diff[j].first << " " << diff[j].second<<endl;
		  }
	      cout << "after set_polars" << endl
	;
   cout<<"at end of constructor" << " hexes counted " << totalHex << " n " << n << endl;
  } //end of constructor

  //function, gets distance between points now supseded by morphologica function
  //float distanceFrom (const pair<float, float> cartesianPoint) const 

  double getdist(point a, point b) {
    double result;
    result = sqrt((a.xval-b.xval)*(a.xval-b.xval) + (a.yval-b.yval)*(a.yval-b.yval));
    return result;
  }

  //function, finds the indexes from original list in sorted order
  // now using stable_sort to ensure consistency in case of hexes equidistant between
  // centres
template <typename T>
vector<int> sort_indexes(const vector<T> &v) {

  // initialize original index locations, iota assigns integer values from 0 onwards
  // to the vector of ints idx
  vector<int> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  stable_sort(idx.begin(), idx.end(),
       [&v](int i1, int i2) {return v[i1] < v[i2];});

  return idx;
}

     int pair2int (pair<int,int> d, int ibase) {
        int result;
        result = d.first*ibase + d.second;
        return result;
    }

    std::pair<int,int> int2pair(int c, int ibase) {
        std::pair<int, int> result;
        result.first = c / ibase;
        result.second = c%ibase;
        return result;
    }

  //test function to ensure hex numbers are not getting scrambled
  void listHex(){
  cout << "in listHex" << endl;
	for (auto h : Hgrid->hexen) {
	 cout << "h.vi " << h.vi << endl;
   }
}   

// set the radialSegments for all regions   
  void setRadialSegments()
  {
    hexGeometry::point a,b;
    for (unsigned int j=0;j<NUMPOINTS;j++) 
	{
      for (unsigned int i=0;i<regionVertex[j].size();i++)
	  {
	    hexGeometry::lineSegment temp;
	    b.first = this->Hgrid->d_x[regionVertex[j][i]];
	    b.second = this->Hgrid->d_y[regionVertex[j][i]];
	    a.first = this->centres[j].first;
	    a.second = this->centres[j].second;
	    temp = hGeo->createLineSegment(a,b);
	    radialSegments[j].push_back(temp);
	  }
	}  
  }//end of method setRadialSegments

    vector<double> getLaplacian(vector<double> Q, double dx) {
        double overds = 1./(1.5*29.0*29.0*dx*dx);
        vector<double> L(n,0.);
        for(auto h : this->Hgrid->hexen){
         int i = int(h.vi);
            L[i]=(Q[N[i][0]]+Q[N[i][1]]+Q[N[i][2]]+Q[N[i][3]]+Q[N[i][4]]+Q[N[i][5]]-6.*Q[i])*overds;
        }
        return L;
    }
        
		vector<double> chemoTaxis(vector<double> Q, vector<double> P, double dx) {
		vector<double> cT(n,0.);
          double overds = 1./(1.5*29.0*29.0*dx*dx);

        for (auto h : Hgrid->hexen) {
		  unsigned int i = h.vi;
        // finite volume method Lee et al. https://doi.org/10.1080/00207160.2013.864392
	      double dr0Q = (Q[N[i][0]]+Q[i])/2.;
	      double dg0Q = (Q[N[i][1]]+Q[i])/2.;
	      double db0Q = (Q[N[i][2]]+Q[i])/2.;
	      double dr1Q = (Q[N[i][3]]+Q[i])/2.;
	      double dg1Q = (Q[N[i][4]]+Q[i])/2.;
	      double db1Q = (Q[N[i][5]]+Q[i])/2.;
	  //double ncentre = Q[i];

          double dr0P = P[N[i][0]]-P[i];
          double dg0P = P[N[i][1]]-P[i];
          double db0P = P[N[i][2]]-P[i];
	      double dr1P = P[N[i][3]]-P[i];
          double dg1P = P[N[i][4]]-P[i];
          double db1P = P[N[i][5]]-P[i];


	  //finite volume for NdelC, h = s/2
	  cT[i] = (dr0Q*dr0P+dg0Q*dg0P+db0Q*db0P+dr1Q*dr1P+dg1Q*dg1P+db1Q*db1P)*overds;

	  //G[i] = (drN*drC+dgN*dgC+dbN*dbC)/(6.*ds*ds) + NN[i]*lapC[i];

        } //matches for on i
		return cT;
	} //end of function chemoTaxis

  //function to timestep coupled equations
    void step(double dt, double Dn, double Dchi, double Dc) {
      dt = dt * 2.5 / Dn;


       // Set up boundary conditions with ghost points
      //cout << " in time step before ghost points" << endl;
      for(auto h : Hgrid->hexen){
	   // cout << "top of ghost loop hex " << h.vi << " x " << h.x << " y " << h.y << endl;
        if(this->Creg[h.vi] > 0){
          for(int j=0;j<6;j++){
		     int i = int(h.vi);
		 // cout << "in ghost loop j = " << j << " i= " << i << " Nbr " << N[h.vi][j] << endl;
             if(N[h.vi][j] == i) {
       	       NN[N[h.vi][j]] = NN[h.vi];
			//   cout << " NN " << NN[N[h.vi][j]] << " NN central " << NN[h.vi] << endl;
	           CC[N[h.vi][j]] = CC[h.vi];
	      }
	     }
	   }
      }



        double beta = 5.;
        double a = 1., b = 1., mu = 1;
		//cout << " before calls to Laplacian " << endl;
        vector<double> lapN = getLaplacian(NN,ds);
        vector<double> lapC = getLaplacian(CC,ds);
        vector<double> cTaxis = chemoTaxis(NN,CC,ds);

        // step N
        for (auto h : Hgrid->hexen) {
          NN[h.vi]+=dt*( a-b*NN[h.vi] + Dn*lapN[h.vi] - Dchi*cTaxis[h.vi]);
        }

        // step C
        double N2;
        for(auto h : Hgrid->hexen){
		    unsigned int i = h.vi;
            N2 = NN[i]*NN[i];
            CC[i]+=dt*( beta*N2/(1.+N2) - mu*CC[i] + Dc*lapC[i] );
        }
    }//end step

  //function to return area of a region
  double regArea (int regNum) {
    double area = 0;
    for (int i=0;i < (int) this->regionIndex[regNum].size();i++){
	  area += 1.;
	}
    return area;
  } //end of funtion regArea


  //function to return perimeter of a region
  double regPerimeter (int regNum) {
    // cout << "in regPerimeter " << endl;
    double perimeter = 0;
    for (auto h : this->regionIndex[regNum])
      if (Creg[h.vi] > 0)
	perimeter += 1.0;
    return perimeter;
  } //end of function regPephieter


//function to return perimeter of a morphed region
double renewRegPerimeter (int regNum) {
  // cout << "in regPerimeter " << endl;
  double perimeter = 0;
  for (auto h : this->regionIndex[regNum])
	if (h.boundaryHex())
      perimeter += 1.0;
  return perimeter;
} //end of function regPephieter


    // transform vector so its mean is zero
    vector<double> meanzero_vector(vector<double> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
        vector <double> result;
        int size = invector.size();
        //meanzero << "size " << size << endl;
        double sum = 0;
        double absSum = 0;
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

	// transform vector of pairs so its mean is (0,0)
	vector<std::pair<double,double>> meanzero_vector(vector<std::pair<double,double>> invector) {
		double sum1, sum2 = 0;
		int size = invector.size();
	    vector<std::pair<double,double>> result;
		for (int i=0;i<size;i++) {
		    sum1 += invector[i].first;
		    sum2 += invector[i].second;
		}
		sum1 = sum1 / (1.0*size);
		sum2 = sum2 / (1.0*size);
        for (int i=0;i<size;i++) {
		    std::pair<double,double> tempPair((invector[i].first - sum1),(invector[i].second - sum2));
		    result.push_back(tempPair);
		}
		return result;

    }
    // return the mean of the absolute values of a  vector
        double absmean_vector(vector<double> invector) {
        //ofstream meanzero ("meanzero.txt",ios::app);
	      double sum = 0;
          int size = invector.size();
        //meanzero << "size " << size << endl;
          for (int i=0; i <size; i++) {
            sum += fabs(invector[i]);
        }
        sum = sum/(1.0*size);
        return sum;
    } 



    double maxVal( vector<double> invector) {
            double result = -1.0e7;
            for (unsigned int i = 0; i < invector.size(); i++) {
                    if (invector[i] > result)
                            result = invector[i];
            }
            return result;
    }


    // to find the minimum value of a vector
    double minVal( vector<double> invector) {
            double result = 1.0e7;
            for (unsigned int i = 0; i < invector.size(); i++) {
                    if (invector[i] < result)
                            result = invector[i];
            }
            return result;
    }
  //function to return fraction of area with NN positive
  double regNNfrac (int regNum) {
    double area = 0;
    double positive_area = 0;
    //int size = this->regionIndex[regNum].size();
    vector<double> normalNN;
    // to normalise the NN field
     for (auto h : this->regionIndex[regNum]){
          normalNN.push_back(this->NN[h.vi]);
      }
      normalNN = this->meanzero_vector(normalNN);
      for (int i=0;i < (int) this->regionIndex[regNum].size();i++){
	  area += 1.;
	  if ((normalNN[i]) > 0) {
		  positive_area += 1.0;
	  }
    }
    return positive_area / area;
  } //end of function regNNfrac

	int printRegion(int regNum) {
	  int result = 0;
	  cout << " in printRegion " << regNum << endl;
        for (auto h : this->regionIndex[regNum]) {
			cout << h.vi<< " ";
			result++;
			}
		return result;
	}

// function to give r and theta relative to region centre
    pair <double,double> set_polars(int regNum){
        pair <double, double> result;
		result.first = 0.0;
		result.second = 0.0;
        double xcentre = this->centres[regNum].first;
        double ycentre = this->centres[regNum].second;
        double xav=0;
        double yav = 0;
        int hexcount = 0;
        cout << "in set polars region " << regNum << " size " << this->regionIndex[regNum].size() << endl;
		/*
  cout << "in listHex" << endl;
	for (auto h : Hgrid->hexen) {
	 cout << "h.vi " << h.vi << endl;
   }
   */
        for (auto h : this->regionIndex[regNum]) {
            hexcount++;
		//	cout << " in set centres hex " << this->regionIndex[regNum][i]<< " i " << i << endl;
            xav += this->Hgrid->d_x[h.vi];
            yav += this->Hgrid->d_y[h.vi];
        }
		//cout << "after set centres " << endl;
		if (hexcount != 0) {
        xav = xav / (hexcount*1.0);
        yav = yav / (hexcount*1.0);
		}
		else {
		  //cout << " in set_polars no hexes in region "<<endl;
		  }
//set the phi values for each hex, this time relative to the region centre
        for (auto&  h : this->regionIndex[regNum]) {
            int index = h.vi;
			double angle;
			//cout <<"in set polars index " << index << " i " << h.vi <<endl;
			cout << "d_x " << this->Hgrid->d_x[index] << " d_y " << this->Hgrid->d_y[index] <<endl;
            h.r = sqrt((this->Hgrid->d_x[index]-xcentre)*(this->Hgrid->d_x[index]-xcentre) 
			+ (this->Hgrid->d_y[index]-ycentre)*(this->Hgrid->d_y[index]-ycentre));
            if ((this->Hgrid->d_y[index] -ycentre) >= 0) {
              angle =  + atan2((this->Hgrid->d_y[index]-ycentre), (this->Hgrid->d_x[index]-xcentre));
			  h.phi = angle;
			  cout<< "region" << regNum << " h.phi "  << h.phi<<  " index " << h.vi << endl;
			  }
            else {
              angle =  2*PI + atan2((this->Hgrid->d_y[index]-ycentre), (this->Hgrid->d_x[index]-xcentre));
              h.phi = angle;
			  cout<< "region " << regNum << " h.phi " << h.phi<<  " index " << h.vi << endl;
			  }
        }
        result.first = xav - xcentre; //diff between seed point and barycentre
        result.second = yav - ycentre ;
		cout << " result.first " << result.first << " result.second " << result.second <<endl;
        return result;
    } //end of function set_polars

    void shift_polars (vector<double> hexVector, double angle) {
        for (auto i = hexVector.begin(); i != hexVector.end();++i)
            hexVector[*i] += angle;
    }


    //function to find all the edges and vertices of the internal boundary
     vector <std::pair<double,double>> dissectBoundary(string logpath) {
         ofstream hfile ( logpath + "/dissectDebug.out" );
         ofstream ifile ( logpath + "/regionList.out" );
         ofstream kfile ( logpath + "/edgesList.out" );
         ofstream lfile ( logpath + "/verticesList.out" );
		 ofstream ufile ( logpath + "/keysList.out" );
         hfile<<"just in dissectBoundary"<<endl;
         vector<std::pair<double,double>> result; 
         vector<int>  regionBoundary; //holding array for the hexes in each region boundary
		 this->edgeIndex.resize(0); //reset edgeIndex
	     int sideCount = 0;
         for (int i=0; i < NUMPOINTS; i++) {
             cout << " region index " << i << " size "  <<this->regionIndex[i].size() << endl;
		     vector<double> psi;
		     psi.resize(0);
             regionBoundary.resize(0);
// fill the regionBoundary
             for (auto h : this->regionIndex[i]) {
		         double angle;
// cout << "boundary i = "<< this->regionIndex[i][j] << " Creg = "<< this->Creg[this->regionIndex[i][j]]<< " j " << j << endl;
                 if (this->Creg[h.vi] >0){
// cout << "boundary i = "<< this->regionIndex[i][j] << " Creg = "<< this->Creg[this->regionIndex[i][j]]<<endl;
                    regionBoundary.push_back(h.vi);
			        angle = h.phi;
			        cout<< " getPhi test " << angle <<  " index " << h.vi << endl;
			        psi.push_back(angle);
                 }
             } //end of loop on a single region

         cout<<"after filling of regionEdge and regionVertex" <<endl;
         cout<<"regionBoundary.size "<<regionBoundary.size() << endl;


         vector<double> rB; //contains the boundary thetas
         vector<int> irB; //holds the sorted boundary indicies
         for (unsigned int j = 0; j < regionBoundary.size();j++){
             rB.push_back(psi[j]);
//cout << " theta value on boundary " << psi[regionBoundary[j]] << endl;
         } //end of loop to fill the boundary angle vector
         irB = sort_indexes(rB); //irB holds sorted indices after sort on theta
         hfile << " irB size " << irB.size() << " rB size " << rB.size() << endl;
		 for (unsigned int i=0;i<irB.size();i++){
		     hfile << " Creg " << Creg[regionBoundary[irB[i]]] << " index " << regionBoundary[irB[i]] << " theta " << rB[irB[i]] <<endl;
		 } //debugging loop for dissectdebug
         unsigned int irBSize = irB.size();
         if (irBSize == 0) {
             cout << "region i " << region[i][0] << " irB size " << irBSize << endl;
             continue;
          } // catches empty boundaries.

// we now walk round the boundary in theta order, find the first vertex then proceed from there to find vertices in order
         unsigned int idissect = 0; //counts number of boundary hexes processed
         unsigned int Vcount = 0; //count of the vertices
         unsigned int Ecount = 0; //count of the edges
         int newVertex; //integer to enumerate the vertices
         unsigned int offset = 0; //number of boundary hexes before the first vertex
         vector<int> ihE; //contains the sorted indicies of each edge
		 //find the offset to the first vertex
         while (offset < irBSize) {
            //cout << " offset " << offset << endl;
            //cout << "Creg " << Creg[regionBoundary[irB[offset]]] <<endl;
             if (Creg[regionBoundary[irB[offset]]] > 1) 
			    {
                     Vcount++; //its a vertex
// cout << " offset " << offset << " idissect " << idissect<<endl;
				     newVertex = regionBoundary[irB[offset]];
				     idissect++;
				     break;
                } 
			    offset++;
		 }
         cout<<"after offset loop" << " offset " << offset  << " idissect " << idissect << endl;
         while ((idissect < irBSize)) {
              Ecount = 0;
              ihE.resize(0);
// while loop to catch the nasty case of artificial adjacent vertices, we only count the last.
// this is a result of our proceeding via hex body rather than hex edge.
              while (Creg[regionBoundary[irB[(idissect + offset) % irBSize]]] > 1) {
		   	      newVertex = regionBoundary[irB[(idissect+offset)%irBSize]];
                  cout << "in vertex loop" << " idissect " << idissect << " Creg "
                  << Creg[regionBoundary[irB[(idissect + offset) % irBSize]]] << endl;
			      idissect++;
                  Vcount++;
			  } // cout << "vertexloop" <<endl;
//walk along the edge until the next vertex
              cout << "Creg " << Creg[regionBoundary[irB[(idissect + offset)%irBSize]]] << " boundary " << regionBoundary[irB[(idissect + offset)%irBSize]] << " newVertex Creg " << Creg[newVertex] << endl;
		      regionVertex[i].push_back(newVertex);
              while ((this->Creg[regionBoundary[irB[(idissect + offset)%irBSize]]] == 1) && (idissect < irBSize)) {
                   ihE.push_back(regionBoundary[irB[(idissect + offset)%irBSize]]);
                   Ecount++;
				   idissect++;
                 }
              ihE.insert(ihE.begin(),newVertex);
              Ecount++;
              cout << "after edge loop Ecount " << Ecount << " Vcount " << Vcount <<endl;
              if (Ecount == 0) {
                  cout<<"Oops - empty edge=========================================="<<endl;
                  continue;
              } //loop to catch empty edge
              cout<<"ihE Size "<<ihE.size() <<endl;
			   
			   /*
			   if (ihE.size() > 1) {
			     cout << " in the ihE loop " << endl;
			     for (unsigned int ii = 0; ii<ihE.size();ii++){
				   int central = region[ihE[ii]][0];
                   for (int ihex = 0; ihex<6; ihex++){ 
				 //  if (hexRegionList[ihE[ii]][ihex] != central) {
                     hfile<<" central " <<region[ihE[ii]][0]<< " nbr region " << hexRegionList[ihE[ii]][ihex]<< " hex " << ihE[ii] << " Creg " << Creg[ihE[ii]] << endl;
					 //central = hexRegionList[ihE[ii]][ihex]; //}
				   }
				   hfile << "$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$4" << endl;
			     }
			   }
			   else {
			     cout << " not in the ihE loop " << endl;
				 continue;
				 }
			   */	 
              cout <<"Vcount "<< Vcount << " Ecount "<< Ecount << endl;
              int regMiddle = region[ihE.rbegin()[1]][0]; // find the first hex in the edge after the vertex edge
              int edgeOuter = -2;
              for (int ihex = 0; ihex<6; ihex++){ //find the first region not the same as the central, since Creg = 1 there can only be one such region.
                  if (regMiddle != hexRegionList[ihE.rbegin()[1]][ihex]){
                      edgeOuter = hexRegionList[ihE.rbegin()[1]][ihex];
                      break;
                  }
              }
              hfile<<"after edgeOuter assignment"<<endl;
              hfile<<"edgeOuter "<<edgeOuter<< " edgeInner " << i << endl;
			  if (edgeOuter > -1) { // edgouter = -1 means that the edge is on the outside of the computational region
                  std::pair <int,int> keypair(i,edgeOuter);
                  int keyint = this->pair2int(keypair,this->base);
                  std::pair <int, vector<int>> p1(keyint,ihE);
                  hfile <<"after pair set keyint " << keyint << " edgeOuter " << edgeOuter << endl;
                  this->edges.insert(p1);
				  if (ihE.size() > 0) {
				  this->edgeIndex.push_back(keyint);
				  sideCount++;
				  }
                  hfile << "after edges insert" << endl;
                  hfile <<"=================================================="<<endl;
                  this->regionList[i].push_back(edgeOuter);
			  }
              else {
                   cout << "edgeouter "<< edgeOuter << endl;
                   this->regionList[i].push_back(edgeOuter);
              } 
          } // end of idissect loop
		  cout << " after idissect loop region " << i << " regionList size " << regionList[i].size()<<  endl;
		  if (regionList[i].size() != Vcount) {
		      cout << "AyAyAy duplicate region" << endl;
		  }
		  if (regionList[i].size() == 0) { 
			  continue;
		  }
//write out to  regionList file the neighbouring regions to this one
		  ifile << "number of nbrs for region " << i << " is " << regionList[i].size() << endl;
          for (unsigned int iregion = 0; iregion < regionList[i].size(); iregion++) {
               ifile << " r " << i << " rNbr " << regionList[i][iregion];
               cout << " r " << i << " rNbr " << regionList[i][iregion];
               ifile << endl;
               ifile << "---------------------------------------"<< endl;
          }
		  cout <<endl;
// write to vertexlist file, vertex x and y coordinates
		  lfile << "number of vertices for region " << i << " is " << regionVertex[i].size() << endl;
          for (unsigned int iregion = 0; iregion < regionVertex[i].size(); iregion++){
               lfile << " r " << i << " vertex " << regionVertex[i][iregion] << " x " << Hgrid->d_x[regionVertex[i][iregion]] << " y " << Hgrid->d_y[regionVertex[i][iregion]] << endl;
               cout << " r " << i << " vertex " << regionVertex[i][iregion];
               lfile << endl;
               lfile << "---------------------------------------"<< endl;
          }
		  cout << " end of iteration of region " << i << endl;

        } //end of loop on regions
		cout << " after loop on regions edgeIndex size " << edgeIndex.size() << " edges size " << edges.size() << " sideCount " << sideCount << endl;
        int difference = 0;
// loops to fill the edges map structure		
        int countIndex = 0;
		int halfway = int(edgeIndex.size()) / 2;
        for (int irlook = 0; irlook < NUMPOINTS; irlook++) {
            for (int jrlook = 0; jrlook < irlook; jrlook++) {
                std::pair <int,int> klook(irlook,jrlook);
                std::pair <int, int> koolk(jrlook,irlook);
                int k = pair2int(klook,this->base);
				int k1 = pair2int(koolk,this->base);
                if (edges.count(k) == 0) // no entry in the container for K
                    continue;
                else {
                    int sizeij = edges[k].size();
                    int sizeji = edges[k1].size();
                    difference = (sizeij - sizeji)*(sizeij - sizeji);
                    if ((sizeij*sizeji != 0) && (difference < 1000)) { // just to check that the edge size hasnt run away
                        kfile << irlook << " " << jrlook << " sizeij " << sizeij << " "<< jrlook << " " << irlook << " sizeji " << sizeji << endl;
				        hfile << " k = " << k << " k1 = " << k1 << " countIndex = " << countIndex << " edge Indexij = " << edgeIndex[countIndex] << " edgeIndexji " <<  edgeIndex[countIndex + halfway] <<  endl;
		                countIndex++;
					}
					else{
					    hfile << " size was zero for k " << k << " and k1 " << k1 << " countIndex " << countIndex << " edge Indexij = " << edgeIndex[countIndex] << " edgeIndexji " <<  edgeIndex[countIndex + halfway]  << endl;
						countIndex++;
					}
                }
            } // end of inner loop over edges
        } //end of outer loop over edges
		for (unsigned int i = 0; i<edgeIndex.size(); i++) {
		    ufile << " for number in edgeIndex " << i << " value is " << edgeIndex[i] << endl;
		}
        return result;
    }//end of function dissectBoundary

    // function to correlate matching edges
    double correlate_edges(string logpath)  
	{
	    double result = 0;     
        ofstream edgefile(logpath + "/edgeCorrelations.txt",ios::app);
        ofstream edgerest(logpath + "/edgeRest.txt",ios::app);
		ofstream correl(logpath + "/correlate.data");
        vector<double> tempvect1;
        vector<double> tempvect2;
        vector<double> tempvect3;
        edgefile << " In correlate_edges " << endl;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        int countResult = 0;
        for (int i = 0; i <NUMPOINTS; i++) 
		{
		    vector<double> normalNN;
		    double NNmean;
		    int size = (int) this->regionIndex[i].size();
		    normalNN.resize(size,0.0);
		    // create a vector of the NN values in the region
		    for (auto h : this->regionIndex[i])
			{
			    //edgefile << "create normalNN j " << j <<endl;
		       	normalNN.push_back(NN[h.vi]);
		    }
	       	NNmean = this->absmean_vector(normalNN);
		    edgefile << " mean of NN in region " << i << "is " <<NNmean << endl;
            edgefile << " i iteration " << i << endl;
            for (auto j = this->regionList[i].begin(); j != this->regionList[i].end();j++) 
			{
                //edgefile << " j iteration " << *j << endl;
                tempvect1.resize(0);
                tempvect2.resize(0);
                tempvect3.resize(0);
                std::pair<int,int> edgePair1(i,*j);
                int edgeIndex1 = this->pair2int(edgePair1,this->base);
                std::pair<int,int>  edgePair2(*j,i);
                int edgeIndex2 = this->pair2int(edgePair2,this->base);
                int count1 = 0;
                int count2 = 0;
                for (auto itr = this->edges[edgeIndex1].begin(); itr != this->edges[edgeIndex1].end();itr++)
				{
                    tempvect1.push_back(this->NN[*itr] - NNmean);
                    count1++;
                }
                // edgefile << " after filling tempvector1 " << endl;
                for (auto itr = this->edges[edgeIndex2].begin(); itr != this->edges[edgeIndex2].end();itr++) 
				{
                     tempvect2.push_back(this->NN[*itr] - NNmean);
                     count2++;
                }
                std::reverse(tempvect2.begin(),tempvect2.end());
                // edgefile << " after filling tempvector2 " << endl;
                //int correlationValue = this->correlate_vector(tempvect1,tempvect2);
                //edgefile << i << " tv1 " << tempvect1.size() << " j " << *j << " tv22  " <<  tempvect2.size() << endl;
                if (tempvect1.size() == tempvect2.size() && tempvect1.size()*tempvect2.size() != 0)
				{
                    double correlationValue = this->correlate_Eqvector(tempvect1,tempvect2);
		            result += fabs(correlationValue);
                    countResult++;
                    edgefile << " i = " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
				    correl << correlationValue << endl;
                    //edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << endl;
                } //end of code if both edges are equal
                else if (tempvect1.size() > tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) 
				{
                    tempvect3 = this->equalize_vector(tempvect2,tempvect1);
                    if (tempvect3.size() == 0)
                        edgefile << " i " << i << " count1 " << tempvect1.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    double correlationValue = this->correlate_Eqvector(tempvect1,tempvect3);
		    result += fabs(correlationValue);
                    countResult++;
                    edgefile << " i " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
				    correl << correlationValue << endl;
                }
                else if (tempvect1.size() < tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) 
				{
                    tempvect3 = this->equalize_vector(tempvect1,tempvect2);
                    if (tempvect3.size() == 0)
                        edgefile << " i " << i << " count2 " << tempvect2.size() << " j " << *j << " tempvect3 is zero  "   << endl;
                    double correlationValue = this->correlate_Eqvector(tempvect3,tempvect2);
		    result += fabs(correlationValue);
                    countResult++;
                    edgefile << " i = " << i << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
				    correl << correlationValue << endl;
                }
            } //end of single edge comparison
        } // end of loop on regions
    edgefile << " countResult "<<countResult<<endl;
	result = result / (countResult * 1.0);
	edgefile.close();
	return result;
    } //end of function correlate_edges

    // method to compare random pairs of edges
	void random_correlate(string logpath, const int max_comp) {
	    ofstream jfile; 
	    ofstream kfile; 
		int max_rand = edges.size();
		jfile.open(logpath + "/random_correlate.txt");
		kfile.open(logpath + "/random_correlate.data");
		vector<double> dinterp, first, second;
		double corr;
		jfile << " max_rand " << max_rand << " max_comp " << max_comp << endl;
		for (int i=0;i<max_comp;i++) {
		   int r1 = rand() % max_rand;
		   int r2 = rand() % max_rand;
		   int rr1 = edgeIndex[r1];
		   int rr2 = edgeIndex[r2]; 
		   int s3 = 0;
		   first.resize(0);
		   second.resize(0);
		   if ((rr1 != 0) && (rr2 != 0))
		   {
		       int s1 = edges[rr1].size();
		       int s2 = edges[rr2].size();
			   jfile << "s1 " << s1 << " s2 " << s2 << endl;
               for (auto itr = this->edges[rr1].begin(); itr != this->edges[rr1].end();itr++) {
                   first.push_back(this->NN[*itr]);
               }
               for (auto itr = this->edges[rr2].begin(); itr != this->edges[rr2].end();itr++) {
                   second.push_back(this->NN[*itr]);
               }
			   first = meanzero_vector(first);
			   second = meanzero_vector(second);
			   if (s1 < s2) {
			       dinterp = equalize_vector(first , second);
				   s3 = dinterp.size();
				   corr = correlate_Eqvector(dinterp, second);
			   }
			   else if (s1 > s2) {
			       dinterp = equalize_vector(second, first);
				   s3 = dinterp.size();
				   corr = correlate_Eqvector(dinterp, first);
			   }
               else {
			       corr = correlate_Eqvector(first, second);
				   s3 = 0;
			   }
               jfile <<  " r1 " << r1 << " rr1 " << rr1 <<" s1 " << s1 << " r2 " << r2 << " rr2 " << rr2 << " s2 " << s2 << " s3 " << s3 << " correlate " << corr << endl;
		       kfile << corr << endl;
		   }
        }
	}

    //function to return the correlation of two vectors
    double correlate_Eqvector(vector<double> vector1, vector<double> vector2) {
        ofstream jfile;
        jfile.open("correlateEqvector.out",ios::app);
        double result;
        jfile << " In correlateEqvector " << endl;
        if (vector1.size() != vector2.size()){
            jfile << "error: vectors must be same length" << endl;
            return -2;
        }
        double vector1Norm = 0;
        double vector2Norm = 0;
        double vector12Product = 0;
        unsigned int vectSize = vector1.size();
        for (unsigned int  i = 0; i != vectSize; i++) {
            vector1Norm +=  vector1[i]*vector1[i];
            vector2Norm += vector2[i]*vector2[i];
            vector12Product += vector1[i]*vector2[i];
        }
            vector1Norm = sqrt(vector1Norm);
            vector2Norm = sqrt(vector2Norm);
            result = vector12Product / (vector1Norm*vector2Norm);
            jfile << "result = " << result << endl;
            return result;

    } //end of function correlateEqvector

    vector <double> equalize_vector (vector<double> svector, vector<double> lvector) {
        ofstream kfile;
        kfile.open("eqVector.out",ios::app);
        vector <double> result;
        int sSize, lSize = 0;
        double sStep, lStep = 0;
        // double delta = 0.000001;
        result.resize(0);
        sSize = svector.size();
        lSize = lvector.size();
        sStep = 1.0 / (1.0 * (sSize-1));
        lStep = 1.0 / (1.0 * (lSize-1));
        double start = 0;
        double finish = 0;
        double value = 0;
        int count = 1;
		double delta = 0.00000001;
        result.push_back(svector[0]);
        for (int i = 0; i < sSize-1; i++) { // walk along the short vector
            start = i*sStep;
            finish = (i+1)*sStep;
            for (int j = 0; j < lSize; j++) { //walk along the long
                if ((j*lStep > start) && (j*lStep) <= finish+delta) {
                    //value = svector[i]*(j*lStep - start) + svector[i+1]*((j+1)*(finish - j*lStep));
                    value = (svector[i]*(j*lStep - start) + svector[i+1]*(finish - j*lStep))/sStep;
                    result.push_back(value);
                    count++;
                }
            }
        }
        if (count != lSize){
           kfile << " count " << count << " lSize " << lSize << " sSize " << sSize << " count " << count <<  " not filled" << endl;
            //  result.resize(0);
            return result;
        }
        else {
             kfile << " count " << count << " l size " << lSize << " sSize " << sSize << " count " << count << " filled" << endl;
            // result.push_back(svector[sSize - 1]);
            return result;
        }
    } //end of function equalize_vector

// determines whether a line segment intersects a hex
	 bool lineSegIntersectHex(hexGeometry::lineSegment s, morph::Hex h) {
	   bool result;
	   bool inside = false;
	   double d = h.getD();
	   if (d != this->ds)
	   {
	     cout << "oops d " << d << "this->ds " <<this->ds<<endl;
	   }	 
	   hexGeometry::point hexCentre;
	   hexCentre.first = h.x;
	   hexCentre.second = h.y;
	   hexGeometry::segDist sDist = this->hGeo->point2Seg(hexCentre,s);
	   inside = sDist.inside;
	   result = (inside  && (sDist.dist <= d/2));
	   if (result) {
       cout << "in lineSegIntersec point " << h.x << " , " << h.y << " h.vi " << h.vi << " distance " << sDist.dist << endl;
	   }
	   return result;
	   }

// determines whether a line segment intersects a hex referred as its index
	 bool lineSegIntersectHex(hexGeometry::lineSegment seg, int index) {
	   bool result;
	   double d = this->ds;
	   hexGeometry::point hexCentre;
	   hexCentre.first = this->Hgrid->d_x[index];
	   hexCentre.second = this->Hgrid->d_y[index];
	   hexGeometry::segDist sDist = this->hGeo->point2Seg(hexCentre,seg);
	   result = (sDist.inside && (sDist.dist <= d/2));
	   /*
	   if (sDist.dist <= d/2)
	   {
	   cout <<" In lineSegIntersectHex dist " << sDist.dist << " d " << d << endl;
	   00}
	   */
	   return result;
	   }
      // method to populate vector of boundary curves
	  void populateBoundVector(bool first) {
        for (int i=0;i<NUMPOINTS;i++) 
		{
		  this->curvedBoundary[i] = this->roundBoundary(i,first);
		}
	  }


// method to round the corners of a region
      morph::BezCurvePath<float> roundBoundary (int regNum, bool first) 
	  {
      morph::BezCurvePath<float> bound;
     // vector<pair<double,double>> vCoords;
     // vector<pair<double,double>> mCoords;
	  int size = this->regionVertex[regNum].size();
	  this->edges.clear();
	  this->vCoords[regNum].resize(size);
	  this->mCoords[regNum].resize(size);
	  //iterate over polygon vertices
	  cout << " Vcoords region " << regNum << endl;
	  if (first) 
	  {
	    for  (int i = 0; i < size;i++)
	    {
	    //make pairs for polygon vertices and their midpoint
	      int vaIndex;
          vaIndex = this->regionVertex[regNum][i];
	    //  vbIndex = this->regionVertex[regNum][(i+1)%size];
		  double firstx = this->Hgrid->d_x[vaIndex];
		  double firsty = this->Hgrid->d_y[vaIndex];
		//  double secondx = this->Hgrid->Hgrid->d_x[vbIndex];
		//  double secondy = this->Hgrid->Hgrid->d_y[vbIndex];
		  pair<double, double> holdPair(firstx,firsty);
	    //  pair2<double, double> holdPair(secondx,secondy);
	      this->vCoords[regNum][i] = holdPair; 
		  //vCoords[i].second = this->Hgrid->d_y[vaIndex];
	     // vCoords[(i+1)%size].first = this->Hgrid->d_x[vbIndex]; 
	     // vCoords[(i+1)%size].second = this->Hgrid->d_y[vbIndex];
		} //end of filling Vcoords
	  } //end region used on first invocation
	  // code for all iterations to round corners
	  for (int i=0;i<size;i++)
	  {  
		this->mCoords[regNum][i].first = (this->vCoords[regNum][i].first + this->vCoords[regNum][(i+1)%size].first)/2.0;
		this->mCoords[regNum][i].second = (this->vCoords[regNum][i].second + this->vCoords[regNum][(i+1)%size].second)/2.0;
	    cout <<" v.x " << this->vCoords[regNum][i].first << " v.y " << this->vCoords[regNum][i].second << " m.x " << this->mCoords[regNum][i].first << "   m.y " << this->mCoords[regNum][i].second << endl;
	  }
	  cout<< endl;
	  //now create the BezCurvePaths
	  for  (int i = 0; i < size;i++)
	  {
	    std::pair<double, double> ma, mb, va;
		ma = this->mCoords[regNum][((i-1)+size)%size];
		mb = this->mCoords[regNum][i];
		va = this->vCoords[regNum][i];
	    morph::BezCurve<float> bc(ma,mb,va);
	    bound.addCurve(bc);
	  }
	  //now update vCoords with mCoords 
	  this->vCoords[regNum] = this->mCoords[regNum];


	  return bound;
	}//end of roundBoundary method
	  
// returns the shortest distace from the seed point to the region boundary
// use of Creg makes it only relevant to unmorphed code
	  double min_radius(int regNum) {
        point  barycentre;
        point boundHex;
        // morphing needs to work with centres, not barycentres
        barycentre.xval =  this->centres[regNum].first;
        barycentre.yval =  this->centres[regNum].second;
        double minradius = 100000.0;
        for (auto h : this->regionIndex[regNum]) {
           // if (Creg[h.vi] > 0) {
		   if (h.boundaryHex() == true)
		   {
                boundHex.xval = Hgrid->d_x[h.vi],
                boundHex.yval = Hgrid->d_y[h.vi];
                double boundDist = getdist(boundHex,barycentre);

                if (boundDist < minradius)
                    minradius = boundDist;
            }
        }
        return minradius;
    }

     double max_radius(int regNum) {
        point  barycentre;
        point boundHex;
		//morphing needs to work with centres, not barycentre
        barycentre.xval = this->centres[regNum].first;
        barycentre.yval = this->centres[regNum].second;
        double maxradius = -100000.0;
        for (auto h : this->regionIndex[regNum]) {
                boundHex.xval = Hgrid->d_x[h.vi],
                boundHex.yval = Hgrid->d_y[h.vi];
                double boundDist = getdist(boundHex,barycentre);

                if (boundDist > maxradius)
                    maxradius = boundDist;
        }
        return maxradius;
    }


    //sectorize over radius
    vector <double> sectorize_reg_radius (int regNum, int numSectors, int beginAngle, int endAngle) {
        ofstream dfile ("logs/sectorRadius.txt",ios::app);
        vector <double>  radiusNN;
        vector <double> normalNN;
	    vector <int> radiusCount;
        radiusCount.resize(numSectors,0); 
	    radiusNN.resize(numSectors,0);
        double startradius, finishradius, radiusInc; //sector radii
        double minradius = min_radius(regNum);
// double maxradius = max_radius(regNum);
        dfile << "region " << regNum << " minradius used " << minradius <<endl;
        radiusInc = minradius /(1.0*numSectors);
        double startAngle, finishAngle, angleInc; //sector angles
        angleInc = 2*PI/(1.*numSectors);
        startAngle = (beginAngle%numSectors)*angleInc;
        finishAngle = (endAngle%numSectors)*angleInc;
 //int size = (int) this->regionIndex[regNum].size();
 // to normalise the NN field
        for (auto h : this->regionIndex[regNum]){
           normalNN.push_back(this->NN[h.vi]);
        }
        normalNN = meanzero_vector(normalNN); 
	    dfile << "after normalise "<<endl;
//for (int i=0;i<size;i++)
// dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

        for (int k=0;k<numSectors;k++) {
        startradius = (k*radiusInc);
        finishradius = (k+1)*radiusInc;
        int count = 0;
        for (auto h : this->regionIndex[regNum]) {
		    if (h.phi >= startAngle && h.phi < finishAngle) {
	            if (h.r >= startradius && h.r <  finishradius) {
	                radiusCount[k]++;
                    radiusNN[k] += normalNN[count];
	            } //end of if on radius
            } //end of if on angleSector
	         count++;
		} //end of loop over the hexes in the region


        if (radiusCount[k] == 0)
    	    radiusNN[k] = 0.0;
        dfile << " region " << regNum << " startradius "<<startradius<<"  finishradius "<<finishradius<< " radiusNN " << radiusNN[k] << endl;

        }//end loop over regions

	    dfile << endl;
        return radiusNN;

    } //end of function sectorize_radius

     //sectorize over radius
    //adapted for digital output
    vector <int> sectorize_reg_Dradius (int regNum, int numSectors, int beginAngle, int endAngle) {
       ofstream dfile ( "logs/sectorRadius.txt",ios::app);
       vector <int>  radiusNN;
       radiusNN.resize(numSectors,0);
       vector <int> radiusCount;
       vector <double> radiusHold;
       vector <double> normalNN;
       radiusCount.resize(numSectors,0);
       radiusHold.resize(numSectors,0);
       double startradius, finishradius, radiusInc; //sector radii
       double maxradius = max_radius(regNum);
       double minradius = min_radius(regNum);
       dfile << "region " << regNum << " maxradius used " << maxradius << " minradius used " << minradius <<endl;
       radiusInc = minradius /(1.0*numSectors);
       double startAngle, finishAngle, angleInc; //sector angles
       angleInc = 2*PI/(1.*numSectors);
       startAngle = (beginAngle%numSectors)*angleInc;
       finishAngle = (endAngle%numSectors)*angleInc;
 //int size = (int) this->regionIndex[regNum].size();
       for (auto h : this->regionIndex[regNum]){
          normalNN.push_back(NN[h.vi]);
       }
// to normalise the NN field
       normalNN = meanzero_vector(normalNN);
       double epsilon = 0.0001*(this->maxVal(normalNN) - this->minVal(normalNN));
 //for (int i=0;i<size;i++)
 // dfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

      for (int k=0;k<numSectors;k++) {
         startradius = (k*radiusInc);
         finishradius = (k+1)*radiusInc;
	     int count = 0;
         for (auto h : this->regionIndex[regNum]) {
            if (h.phi >= startAngle && h.phi < finishAngle) {
                if (h.r >= startradius && h.r < finishradius) {
                    radiusCount[k]++;
//radiusCC[k] += this->CC[this->regionIndex[regNum][i]];
                   radiusHold[k] += normalNN[count];
                } //end of if on radius
            } //end of if on angleSector
            count++;
         } //end of loop over hexes in and individual region
      }//end of loop over all regions
	
	  dfile << "after creation of sectorized field region Dregion " << regNum <<  endl;
    
      for (int k=0;k<numSectors;k++){
	     startradius = (k*radiusInc);
	     finishradius = (k+1)*radiusInc;
	     if (radiusCount[k] == 0) {
		    radiusNN[k] = 2;
		    continue;
	     }
//radiusHold[k]  = radiusHold[k]  / (1.*radiusCount[k]);
         if (radiusHold[k] > epsilon)
             radiusNN[k] = 1;
         else if (radiusHold[k] < epsilon)
             radiusNN[k] = -1;
         else
             radiusNN[k] = 0;

         dfile << " region " << regNum <<" startradius "<<startradius<<"  finishradius "<<finishradius<< " DradiusNN " << radiusNN[k] << endl;
      }//end loop over sectors
	  dfile << endl;
      return radiusNN;

   } //end of function sectorize_radius



 //function to count the hexes in sectors of a region via angular sectors
    vector <double> sectorize_reg_angle (int regNum, int numSectors, int beginradius, int endradius) {
    //std::pair<double,double> diff; //difference between seed point and CoG of region
        ofstream cfile ("logs/sectorAngle.txt",ios::app);
        vector <double> angleNN; //average value of CC in each sector
        vector <double> normalNN;
        vector <int> angleCount; //number of hexes in each sector
        angleNN.resize(numSectors,0);
        angleCount.resize(numSectors,0);
        double startAngle, endAngle, angleInc; //sector angles
        double startradius, finishradius,radiusInc;
        double minradius = min_radius(regNum);
// double minradius = min_radius(regNum);
        radiusInc = minradius/ (1.0*numSectors);
        startradius = beginradius*radiusInc;
        finishradius = endradius*radiusInc;
        angleInc = 2*PI/(1.*numSectors);
// to normalise the NN field
//int size = (int) this->regionIndex[regNum].size();
        for (auto h : this->regionIndex[regNum]){
           normalNN.push_back(this->NN[h.vi]);
        }
        normalNN = meanzero_vector(normalNN);
//for (int i=0;i<size;i++)
//   cfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

        for (int k=0;k<numSectors;k++) {
            startAngle = k*angleInc;
            endAngle = (k+1)*angleInc;
            if ((k+1) == numSectors)
               endAngle = 2*PI;

            int count = 0;
            for (auto h : this->regionIndex[regNum]) {
                if ( h.r  >= startradius && h.r < finishradius) {
                    if (h.phi >= startAngle && h.phi < endAngle) {
                        angleCount[k]++;
//angle[k] += this->[this->regionIndex[regNum][i]];
                        angleNN[k] += normalNN[count];
                    }//end if on angle
//cfile << setw(5) << angleVal[i]  <<"  ";
                }//end if on radius
		        count++;
            }//end loop over all hexes in a region
        }//end loop on all sectors
	    cfile << "after creation of sectorized field region angle " << regNum << " number of hexes " << count <<  endl;

        angleNN = meanzero_vector(angleNN);
        for (int k=0;k<numSectors;k++){ //calculate the average angle in the sector
	        startAngle = k*angleInc;
	        endAngle = k*angleInc;
	        if (angleCount[k] != 0)
		        angleNN[k] = angleNN[k]/(1.*angleCount[k]);
	        else
		        angleNN[k] = -999.999;
//write out values
		    cfile << " region " << regNum <<" startAngle "<< startAngle << "  endAngle "<< endAngle << " angleNN " << angleNN[k] << endl;
        }//end loop on sectors

        cfile << endl;
        return angleNN;

    } //end of function sectorize_region

 //function to count the hexes in sectors of a region via angular sectors
    vector <int> sectorize_reg_Dangle (int regNum, int numSectors, int beginradius, int endradius) {
        ofstream cfile ("logs/sectorAngle.txt",ios::app);
 //std::pair<double,double> diff; //difference between seed point and CoG of region
        vector <int> angleNN; //digitized value of NN in each sector
        vector <double> angleHold;
        vector <double> normalNN;
        vector <int> angleCount; //number of hexes in each sector
        angleNN.resize(numSectors,0);
        angleHold.resize(numSectors,0);
        angleCount.resize(numSectors,0);
        double startAngle, endAngle, angleInc; //sector angles
        double startradius, finishradius,radiusInc;
 //double maxradius = max_radius(regNum);
        double minradius = min_radius(regNum);
        radiusInc = minradius/ (1.0*numSectors);
        startradius = beginradius*radiusInc;
        finishradius = endradius*radiusInc;
        angleInc = 2*PI/(1.*numSectors);
// to normalise the NN field
//int size = (int) this->regionIndex[regNum].size();
    	int NNcount = 0;
        for (auto h : this->regionIndex[regNum]){
            normalNN.push_back(this->NN[h.vi]);
		    cfile << " h.vi " << h.vi << " NNcount " << " h.phi " << h.phi << " h.r " << h.r << NNcount <<  endl;
		    NNcount++;
        }
        normalNN = meanzero_vector(normalNN);
        double epsilon = 0.0001*(this->maxVal(normalNN) - this->minVal(normalNN));
	    cfile << " after normalisation NN field in Dangle region " << regNum <<  " NNcount " << NNcount << endl;
// for (int i=0;i<size;i++)
//    cfile << " i " << i << " normalNN[i] " << normalNN[i] << endl;

        for (int k=0;k<numSectors;k++) {
            startAngle = k*angleInc;
            endAngle = (k+1)*angleInc;
            if ((k+1) == numSectors)
            endAngle = 2*PI;
            cfile << " start of numSectors loop " << k << endl;		 
            int count = 0;
            for (auto &h : this->regionIndex[regNum]) {
                cfile << " start of region index  loop " << count <<  " hex " << h.vi << " h.phi " << h.phi << " h.r " << h.r << endl;		 
                if (h.r >= startradius && h.r < finishradius) {
                    if (h.phi >=startAngle && h.phi < endAngle) {
//angleVal.push_back(h.phi);
                    angleCount[k]++;
//angle[k] += this->[this->regionIndex[regNum][i]];
                    angleHold[k] += normalNN[count];
                    }//end if on angle
//cfile << setw(5) << angleVal[angleCount[k]]  <<" region  " << regNum << endl;
                }//end if on radius
	            count++;
            } //end if over hexes in a region
        }//end if on over all sectors
	    cfile << "after creation of sectorized field region Dangle " << regNum <<  endl;

        angleHold = meanzero_vector(angleHold);

        for (int k=0;k<numSectors;k++){
            startAngle = k*angleInc;
            endAngle = (k+1)*angleInc;
            if (angleCount[k] == 0) {
                angleNN[k] = 2;
                continue;
            }
//angleHold[k] = angleHold[k] / (1.*angleCount[k]);
            if (angleHold[k] > epsilon)
                angleNN[k] = 1;
            else if (angleHold[k] < epsilon)
                angleNN[k] = -1;
            else
                angleNN[k] = 0;

            cfile << " region " << regNum <<" startangle  " << startAngle << "  endAngle  "<< endAngle << " DangleNN " << angleNN[k] << endl;
        } //end loop over sectors
        return angleNN;

    } //end of function sectorize_region digital version

  //method for comparing hexes by angle
  bool hexcompare(morph::Hex h1, morph::Hex h2)
  {
    bool result;
	result = (h1.phi >= h2.phi);
	return result;
  }


  //method to renew a region after rounding
  void renewRegion(int regNum, list<Hex> hexen)
  {
    //regionIndex[regNum].resize(0);
	for (auto& h : hexen) //repopulate regionIndex
	{
	  this->regionIndex[regNum].push_back(h);
    }
  }
   
  //method to renew polars and boundary
  void renewBoundary(int regNum, list<Hex> hexen)
  {
    //pair<double,double> diff;
	for (auto h : hexen) 
	{
      if (h.boundaryHex())
	  {
	    this->regionBound[regNum].push_back(h);
	  }	
    }
	cout << " region " << regNum << " bound size " <<regionBound[regNum].size() << endl;
    //result.first = diff.first + centres[regNum].first;
    //result.second = diff.second + centres[regNum].second;
	return;
  }
   
// now sort the boundary by angle
  void renewDissect(int regNum)
  {   
    vector<int> boundary; //contains the boundary h.vi values
    vector<double> rB; //contains the boundary thetas
    vector<int> irB; //holds the sorted boundary indicies
//	hexGeometry::point a,b;
	int bsize = regionBound[regNum].size();
	int vsize = regionVertex[regNum].size();
	vector<double> vertexAngle;
    
    for (auto h : regionBound[regNum])
	  {
	    double angle = h.phi;
        rB.push_back(angle);
		boundary.push_back(h.vi);
		cout << "region " << regNum <<" theta boundary " << angle << " boundary index " << h.vi << endl;
      }
      irB = sort_indexes(rB); //indices after sort on theta
        cout << " irB size " << irB.size() << " rB size " << rB.size() << endl;
	  sortedBoundary[regNum] = irB;	
	  //print out the line segments
	  for (int i=0;i<vsize;i++){
	    double startx = radialSegments[regNum][i].start.first;
	    double starty = radialSegments[regNum][i].start.second;
	    double endx = radialSegments[regNum][i].end.first;
	    double endy = radialSegments[regNum][i].end.second;
		if (endy >= starty) 
		{
		  vertexAngle.push_back(atan2((endy - starty) , (endx - startx)));
		}
		else 
		{
		  vertexAngle.push_back(2*PI + atan2((endy - starty) , (endx - startx)));
		}
		cout << " start.x " <<radialSegments[regNum][i].start.first  << " start.y " << radialSegments[regNum][i].start.second << " end.x "<< radialSegments[regNum][i].end.first << " end.y " <<radialSegments[regNum][i].end.second<<" angle " << vertexAngle[i] << endl;
      }
	   //write the indices in phi order
	  for (int i = 0; i < bsize; i++) 
	  {
	    cout << " boundHex " << i << " hex " << irB[i] << " angle " << rB[irB[i]] << endl;
	  }
	  int offset = 0;
	  int idissect = 0;
	  vector<vector<int>> ihE;
	  ihE.resize(vsize);
	  while ((rB[irB[offset]] < vertexAngle[0]) && (offset<bsize)) 
	  {
	    cout << " offset inside " << offset << " vertexAngle " << vertexAngle[0] << " hex angle " << rB[irB[offset]]<< " hex " << irB[offset] << endl;
		offset++;
	  } 
	  cout << " offset outside " << offset << " hex angle " << rB[irB[offset]] << " vsize " << vsize << " hex " << irB[offset] << endl;
	  idissect = offset;
	  for (int i=1; i< vsize; i++)
	  {
	    cout << "head of segment loop " << i << " idissect  " << idissect << " vertexAngle " << vertexAngle[i%vsize] << " hex " << irB[offset] << endl;
		  while ((rB[irB[idissect%bsize]] < vertexAngle[i%vsize]) && (idissect <= bsize))
		  {
		    cout << " filling edge loop " << i-1 << " index " << idissect << " angle " << rB[irB[idissect%bsize]] << " hex " << irB[offset%bsize] << endl;
		    ihE[i-1].push_back(irB[idissect%bsize]);
			idissect++;
		  }
		}
		cout << " just before  vsize end loop angle " << vertexAngle[0] + 2*PI << endl; 
		while ((rB[irB[idissect%bsize]] < vertexAngle[0] + 2*PI) && (idissect< bsize+offset-1))
		{
		  cout << " filling end edge loop 3 " << idissect << " angle " << rB[irB[idissect%bsize]] << " hex " << irB[idissect%bsize] <<  endl;
		  ihE[vsize-1].push_back(irB[idissect%bsize]);
		  idissect++;
		}

	  // int index = 0;
	  
	   for (int i=0;i<vsize;i++)
	   {
	     for (unsigned int j=0; j<ihE[i].size(); j++) 
		 {
		  cout << " region " << regNum << " edge " << i << " size " << ihE[i].size() << " index "<< ihE[i][j] << " phi " << rB[ihE[i][j]] << endl;
		//  index++;
		 }
	   }
      
	 for (int iregion = 0;iregion<vsize;iregion++) 
	 {
	   int edgeOuter = regionList[regNum][iregion];
	   cout << "edgeOuter " << edgeOuter << " region " << regNum << endl;
	   if (edgeOuter > -1) 
	   {
	     std::pair<int,int> keypair(regNum,edgeOuter);
		 int keyint = this->pair2int(keypair,this->base);
         std::pair <int, vector<int>> p1(keyint,ihE[iregion]);
       //  cout <<"after pair set"<<endl;
         this->edges.insert(p1);
       //  cout << "after edges insert" << endl;
       //  cout <<"=================================================="<<endl;
	   }
       else 
	   {
         //cout << "edgeouter "<< edgeOuter << endl;
       } 
     } // end filling edges loop

	 //now print out the edges
	 int difference = 0;
	 for (int jrlook = 0; jrlook < regNum; jrlook++) 
	 {
       std::pair <int,int> klook(regNum,jrlook);
       int k = pair2int(klook,this->base);
       if (edges.count(k) == 0)
         continue;
       else 
	   {
         int sizeij = edges[k].size();
         std::pair <int, int> koolk(jrlook,regNum);
         k = pair2int(koolk,this->base);
         int sizeji = this->edges[k].size();
         difference = (sizeij - sizeji)*(sizeij - sizeji);
         if ((sizeij*sizeji != 0) && (difference < 1000))
           cout << regNum << " " << jrlook << " sizeij " << sizeij << " "<< jrlook << " " << regNum << " sizeji " << sizeji << endl;
		 }
	   }
                

   } //end of method renewDissect
   

    // function to renew correlate matching edges
     double renewcorrelate_edges(int regNum, string logpath)  
	  {
	    double result = 0;     
        ofstream edgefile(logpath + "/edgeCorrelations.txt",ios::app);
    	edgefile << " in morphed edge correlation routine "<<endl;
        vector<double> tempvect1;
        vector<double> tempvect2;
        vector<double> tempvect3;
        cout << " In renewcorrelate_edges " << endl;
        // iterate over regions
        //for each region iterate over region edges
        // for each edge pair (i,j), (j,i) call correlate_vectors
        // write the i,j and correlation to a file
        // what happens if there are multiple entries in regionList at the start and the end?
        // there are some regions that have this
        int countResult = 0;
		vector<double> normalNN;
		double NNmean;
		int size = (int) this->regionIndex[regNum].size();
		normalNN.resize(size,0.0);
		// create a vector of the NN values in the region
		for (auto h : this->regionIndex[regNum])
		{
			//edgefile << "create normalNN j " << j <<endl;
		  normalNN.push_back(NN[h.vi]);
		}
	    NNmean = this->absmean_vector(normalNN);
		cout << " mean of NN in region " << regNum << "is " <<NNmean << endl;
        //edgefile << " i iteration " << i << endl;
        for (auto j = this->regionList[regNum].begin(); j != this->regionList[regNum].end();j++) 
		{
          //edgefile << " j iteration " << *j << endl;
          tempvect1.resize(0);
          tempvect2.resize(0);
          tempvect3.resize(0);
          std::pair<int,int> edgePair1(regNum,*j);
          int edgeIndex1 = this->pair2int(edgePair1,this->base);
          std::pair<int,int>  edgePair2(*j,regNum);
          int edgeIndex2 = this->pair2int(edgePair2,this->base);
          int count1 = 0;
          int count2 = 0;
          for (auto itr = this->edges[edgeIndex1].begin(); itr != this->edges[edgeIndex1].end();itr++)
		  {
            tempvect1.push_back(this->NN[*itr] - NNmean);
            count1++;
          }
          // edgefile << " after filling tempvector1 " << endl;
          for (auto itr = this->edges[edgeIndex2].begin(); itr != this->edges[edgeIndex2].end();itr++)
		  {
            tempvect2.push_back(this->NN[*itr] - NNmean);
            count2++;
          }
          std::reverse(tempvect2.begin(),tempvect2.end());
          // edgefile << " after filling tempvector2 " << endl;
          //int correlationValue = this->correlate_vector(tempvect1,tempvect2);
          //edgefile << i << " tv1 " << tempvect1.size() << " j " << *j << " tv22  " <<  tempvect2.size() << endl;
          if (tempvect1.size() == tempvect2.size() && tempvect1.size()*tempvect2.size() != 0)
		  {
            double correlationValue = this->correlate_Eqvector(tempvect1,tempvect2);
		    result += fabs(correlationValue);
            countResult++;
            edgefile << " regNum = " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
            //edgefile << i << " Size1 " << edges[edgeIndex1].size() << " j " << *j << " Size2  " << edges[edgeIndex2].size() << endl;
            } //end of code if both edges are equal
            else if (tempvect1.size() > tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) 
			{
              tempvect3 = this->equalize_vector(tempvect2,tempvect1);
              if (tempvect3.size() == 0)
			  {
                edgefile  << " regNum " << regNum << " count1 " << tempvect1.size() << " count2 " << tempvect2.size() << " j " << *j << " tempvect3 is zero  "   << endl;
			   }
              double correlationValue = this->correlate_Eqvector(tempvect1,tempvect3);
		      result += fabs(correlationValue);
              countResult++;
              edgefile << " regNum " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
            }
            else if (tempvect1.size() < tempvect2.size() && tempvect1.size()*tempvect2.size() != 0) 
			{
              tempvect3 = this->equalize_vector(tempvect1,tempvect2);
              if (tempvect3.size() == 0)
			  {
                edgefile << " regNum is " << regNum << " count2 " << tempvect2.size() << " j " << *j << " tempvect3 is zero  "   << endl;
		      }
              double correlationValue = this->correlate_Eqvector(tempvect3,tempvect2);
		      result += fabs(correlationValue);
              countResult++;
              edgefile << " regNum = " << regNum << " c1 " << count1 << " j " << *j << " c2  " <<  count2 << " cV " << correlationValue << endl;
             }
      } //end of single edge comparison
	  cout << " countResult "<<countResult<<endl;
	  result = result / (countResult * 1.0);
	  edgefile.close();
	  return result;
     } //end of function renewcorrelate_edges

}; // DRegion
