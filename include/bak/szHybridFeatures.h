#ifndef ___SZ_HYBRID_FEATURES_H___
#define ___SZ_HYBRID_FEATURES_H___

#include <vector>
using namespace std;
#include <szMexUtility.h>

inline int
NumHybridFeatures();

vector<CFeature>
computeHybridFeatures(const vector<unsigned short>& A, //image
					  const vector<unsigned char>& B, //nodule segmentation
					  const vector<float>& D, //
					  //int x, int y, int z,
					  const int* dims);


#endif /* ___SZ_HYBRID_FEATURES_H___ */