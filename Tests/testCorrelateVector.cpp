// driver program to test the correlation of two vectors
#include "region.h"
#include "analysis.h"
#include "ksSolver.h"
#include <morph/display.h>
#include <cctype>
#include <locale>
#include <algorithm>
#include <string>
using std::array;
using std::string;
using std::stringstream;
using std::cerr;
using std::endl;
using std::runtime_error;
/*
using morph::HexGrid;
using morph::HdfData;:
using morph::Tools;
using morph::Display
*/
using namespace morph;
using namespace std;
#define N2 100
#define N1 50
#define M 2000

int main (int argc, char **argv)
{
    ofstream afile ("./correlateVector.out");
// seed the random number generator
    unsigned int a_seed = time(0);
	DRegion R(8,"./logs");
	srand(a_seed);
    vector<double> A,A1,B,C;
	A.resize(N1,0.0);
	B.resize(N2,0.0);
	C.resize(M,0.0);
    for (int j=0; j<M; j++) {   
        for (int i=0; i< N1; i++) {
            double choice = morph::Tools::randDouble();
		    if (choice > 0.5)
			{
                A[i] = -morph::Tools::randDouble();
			}	
			else
			{
                A[i]  = morph::Tools::randDouble();
			}
	     }//end of loop over A 
		 for (int i=0; i< N2; i++) {
            double choice = morph::Tools::randDouble();
		    if (choice > 0.5)
			{
                B[i] = -morph::Tools::randDouble();
			}	
			else
			{
                B[i]  = morph::Tools::randDouble();
			}
	    }//end of loop over N
		A1 = R.equalize_vector(A,B);
        C[j] = R.correlate_Eqvector(A1,B);
	} //end of loop over M
	for (int i=0;i<M;i++) {
	   afile << C[i] << endl;
	}   
} //end of main program
   
    
