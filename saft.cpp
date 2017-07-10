#include <iostream>
#include <typeinfo>
#include <string>
#include <armadillo>
#include <time.h>
#include "Headers/myMath.h"
#include "Headers/myArmadillo.h"
#include "Headers/fileoperation.h"

// COMPILE WITH g++ -O3 file.cpp -larmadillo -o saft
// RUN WITH ./saft file.xyz 10 1 100 5.0

using namespace std;
using namespace arma;

int main(int argc, char** argv){

  mat TessCoeff(2601,3);
  TessCoeff.load("DATs/TessarCoeff.dat");

  int lMax = strtol(argv[2], NULL, 10);
  int rMin = strtol(argv[3], NULL, 10);
  int rMax = strtol(argv[4], NULL, 10);
  double sig = atof(argv[5]);
  sig = 1/sig;
  sig = sig*sig;
  
  mat A;
  A.load("TXTs/parameters100.txt");
//	A.load("parameters50.txt");
//   A.load("P200_both.txt");

  mat coord = getPos(argv[1]);
  string* type = getType(argv[1]);

  vec typeVal = sqrt(getTypeVal(type, coord));
//  vec typeVal = sqrt(sqrt(getTypeVal(type, coord)));
//
//  typeVal.print("TypeVals");

  coord = posAve(coord);
  double k;
  double val = 0;
  for(int r=rMin; r <= rMax; r++) { 
     k = 0;
	   	for(int l=0; l <= lMax; l++) {
            for(int m=-l; m <= l; m++){ 

			  	val = tesseralInt2dTableType(l, m, (r + 1)/10.0 , TessCoeff(k,2), coord, typeVal, sig, A) ;
                cout << val ;
                cout << " " ;
			  	k+=1;

	      }
	    }
      cout << endl;
    }

	return 0;
}

