#ifndef ___SZ_EXTRACT_FOREGROUND_H___
#define ___SZ_EXTRACT_FOREGROUND_H___

#include <vector>
using namespace std;

int
doForegroundExtraction(const vector<unsigned short>& A, //input image
					   vector<unsigned char>& L, //output - foreground segmentation 
					   vector<float>& Q, //output - probability image
					   vector<float>& D, //output - distance map image
					   vector<float>& Sp, //output - sphericity image
					   vector<unsigned char>& P, //output - sphericity local maxima
					   vector<int>& vseeds, //output - seed points
					   vector<float>& vsph, //sphericity values at the seed points
					   const vector<int>& vnoduleTypes, //input - nodule types of interest
					   //what follow are various parameters to the RD segmentation
					   int RDIter,
					   float RDlambda,
					   //what follo are various parameters to the seed search
					   int findSeedWidth, 
					   int findSeedCoreWidth,
					   float findSeedThres,
					   //a parameter that enables/disables persistent segmentation
					   bool bPersistent,
					   //dimension of the data
					   int ndim,
					   const int* dims);

#endif /* ___SZ_EXTRACT_FOREGROUND_H___ */