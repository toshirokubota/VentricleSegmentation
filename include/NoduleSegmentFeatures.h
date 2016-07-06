#ifndef ___NODULE_SEGMENT_FEATURES_H___
#define ___NODULE_SEGMENT_FEATURES_H___

#include <vector>
using namespace std;
#include <szFeatureClass.h>

vector<CFeature>
computeFeatuers(const vector<unsigned short>& A, //original data
				const vector<unsigned char>& B,	 //segmentation
				const vector<unsigned char>& B0, //free-trace segmentation
				const vector<unsigned char>& L,  //foreground
				const vector<int>& seed,
				const vector<float>& v,
				int ndimA, 
				const int* dimsA
				);

#endif /* ___NODULE_SEGMENT_FEATURES_H___ */
