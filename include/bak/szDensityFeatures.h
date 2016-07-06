#ifndef ___SZ_DENSITY_FEATURES_H___
#define ___SZ_DENSITY_FEATURES_H___

#include<vector>
using namespace std;

inline int
NumDensityFeatures();

#include <szMexUtility.h>
vector<CFeature>
computeDensityFeatures(const vector<unsigned short>& A, //image
					   const vector<unsigned char>& B, //nodule segmentation
					   int ndim,
					   const int* dims);

#endif /* ___SZ_DENSITY_FEATURES_H___ */