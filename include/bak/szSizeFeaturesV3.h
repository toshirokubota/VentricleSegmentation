#ifndef ___SZ_SIZE_FEATURES_V3_H___
#define ___SZ_SIZE_FEATURES_V3_H___

#include<vector>
using namespace std;

inline int
NumSizeFeatures();

vector<CFeature>
computeSizeFeaturesV3(const vector<unsigned char>& B, //nodule segmentation
					  const vector<float>& v,         //voxel size
					  int ndim,
					  const int* dims);

#endif /* ___SZ_SIZE_FEATURES_V3_H___ */