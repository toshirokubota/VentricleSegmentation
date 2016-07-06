#ifndef ___CAD_NODULE_SEGMENT_RUN_H___
#define ___CAD_NODULE_SEGMENT_RUN_H___

#include <vector>
using namespace std;
#include <szFeatureClass.h>

bool
NoduleSegmentFeatures(vector<CFeature>& vFeatures,
					  const vector<unsigned short>& A,
					  const vector<float>& v,
					  int ndimA,
					  const int* dimsA
					  );

#endif /* ___CAD_NODULE_SEGMENT_RUN_H___ */
