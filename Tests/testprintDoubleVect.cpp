#include <vector>
#include <iostream>
#include <fstream>

using namespace std;

//function to print a vector
void printDoubleVect (std::string path, vector<double> invect) {
    ofstream Vout (path ,ios::app);
    unsigned int size = invect.size();
    vector<double>::iterator ptr;
    Vout << " invect size " << size << endl;
    for (ptr = invect.begin(); ptr < invect.end(); ptr++) {
       Vout << " vector value " << *ptr << "  ";
    }
    Vout << endl << "end of Vector " << endl;
}

 int main() {
    vector<double> v;
    v.resize(0, 0.0);
    for (int i=0; i<10; i++) {
       double temp = i*9.789;
       v.push_back(temp);
    }
    printDoubleVect("./testDoubleVect.txt",v);
    return 0;
}

