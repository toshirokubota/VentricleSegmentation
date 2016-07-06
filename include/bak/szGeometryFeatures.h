#ifndef ___SZ_GEOMETRY_FEATURES_H___
#define ___SZ_GEOMETRY_FEATURES_H___

#include<vector>
using namespace std;

inline int
NumGeometryFeatures();

#include <szMexUtility.h>
vector<CFeature>
computeGeometryFeatures(const vector<unsigned char>& B, //nodule segmentation
						const vector<unsigned char>& B0, //free trace segmentation
						const vector<unsigned char>& L, //foreground segmentation
						const vector<float>& D,         //distance map of B
						const vector<float>& D0,         //distance map of B0
						const vector<float>& Dl,         //distance map of L
						int xseed, int yseed, int zseed, //seed points
						int ndim,
						const int* dims);

vector<CFeature>
computeGeometryFeaturesNew(const vector<unsigned char>& B, //nodule segmentation
						   const vector<unsigned char>& B0, //fee trace result
						   const vector<unsigned char>& L, //foreground including the segmentation
						   const vector<float>& D,
						   const vector<float>& D0,
						   const vector<float>& Dl,
						   int xseed, int yseed, int zseed, //seed points
						   const int* dims);

vector<CFeature>
computeErosionFeatures(const vector<unsigned char>& A, //segmentation
					   const vector<unsigned char>& B0, //free trace segmentation
					   const vector<float>& D0,			//distance map of free trace segmentation
					   const int* dims);

#endif /* ___SZ_GEOMETRY_FEATURES_H___ */