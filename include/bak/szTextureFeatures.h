#ifndef ___SZ_TEXTURE_FEATURES_H___
#define ___SZ_TEXTURE_FEATURES_H___

#include<vector>
using namespace std;

inline int
NumTextureFeatures();

vector<CFeature>
computeTextureFeatures(const vector<unsigned short>& A, //image
					   const vector<unsigned char>& B, //nodule segmentation
					   int ndim,
					   const int* dims);

#endif /* ___SZ_TEXTURE_FEATURES_H___ */