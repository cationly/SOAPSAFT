#include <iostream>
#include <string>
#include <armadillo>
#include <time.h>
#include "myMath.h"
#include "myArmadillo.h"
#include "fileoperation.h"
// COMPILE WITH g++ -O3 file.cpp -larmadillo
using namespace std;
using namespace arma;

int main(int argc, char* argv[]){

mat TessCoeff(2601,3);
TessCoeff.load("TessarCoeff.dat");

	mat A;
	A.load("parameters100.txt");
//	A.load("P200_both.txt");

mat coord1 = getPos(argv[1]);
mat coord2 = getPos(argv[2]);
//mat coord2 = coord1;
  int lMax = strtol(argv[3], NULL, 10);
  double sig = atof(argv[4]);

string* type1 = getType(argv[1]);
//string* type2 = getType("dsgdb9nsd_071702.xyz");
string* type2 = getType(argv[2]);
//type2[0] = "H";

vec typeVal1 = sqrt(sqrt(getTypeVal(type1, coord1)));
vec typeVal2 = sqrt(sqrt(getTypeVal(type2, coord2)));

typeVal1.print("Type");

	coord1 = posAve(coord1);
	coord2 = posAve(coord2);

//  coord2 = rotate3d(coord2,1.0, 1.0,1.0);
//  coord1(51,2) = coord1(51,2) + 0.1;
//    typeVal1(0) = 1.000; 

	double val1 = 0;
	double val2 = 0;
	double accum1 = 0;
	double accum2 = 0;
	double rDiff = 0;
  double radiuses[6] = {0.125,0.25,0.5,1.0,2.0, 4.0};
//  double radiuses[5] = {0.1,0.2,0.5,1.0,2.0};

	for(int r=0; r < 6; r++) { 
		int k = 0;
		for(int l=0; l <= lMax; l = 1 + l) {
			for(int m=-l; m <= l; m++) { 
				val1 = tesseralInt2dTableType(l, m, radiuses[r], TessCoeff(k,2), coord1, typeVal1, sig, A);
				val2 = tesseralInt2dTableType(l, m, radiuses[r], TessCoeff(k,2), coord2, typeVal2, sig, A);
				accum1 += val1*val1;
				accum2 += val2*val2;
				k+=1;

//      cout <<"l: " << l << " r: " << r  << " ms: " << val1 << " " << val2<< endl;

//        cout << l << endl;
//        cout << m << endl;
//        cout << k << endl;
				
			}

      cout <<"l: " << l << " r: " <<radiuses[r]  << " ms: " << sqrt(accum1) << " " << sqrt(accum2)<< endl;


			rDiff = rDiff + abs(accum1 - accum2);
			accum1 = 0;
			accum2 = 0; 
		}

	}
	accum1 = 0; accum2 = 0;
	cout << rDiff << endl;

	return 0;
}
