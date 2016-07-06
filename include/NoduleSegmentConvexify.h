#ifndef ___NODULE_SEGMENT_CONVEXIFY_H___
#define ___NODULE_SEGMENT_CONVEXIFY_H___

#include <vector>
using namespace std;

bool
convexitySegmentation(vector<unsigned char>& B,		  //result
					  const vector<unsigned char>& A, //segmentation
					  const vector<unsigned char>& L, //foreground
					  int ndimA,
					  const int* dimsA);

#endif /* ___NODULE_SEGMENT_CONVEXIFY_H___ */
