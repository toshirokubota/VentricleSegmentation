
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>
using namespace std;

#include <szMiscOperations.h>
#include <szMexUtility.h>	
#include <szMexUtilityTemplate.h>
#include <szConvolutionTemplate.h>
#include <szFeatureClass.h>

inline int
NumTextureFeatures()
{
	const int NumFeatures = 6;
	return NumFeatures;
}

vector<CFeature>
computeTextureFeatures(const vector<unsigned short>& A, //image
					   const vector<unsigned char>& B, //nodule segmentation
					   int ndim,
					   const int* dims) 
{
	vector<CFeature> features(NumTextureFeatures());

	return features;
}
