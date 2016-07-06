#ifndef ___SZ_EXTRACT_NODULE_H____
#define ___SZ_EXTRACT_NODULE_H____

#include <vector>
using namespace std;

int
doNoduleExtraction(const vector<unsigned char>& L, //input - foreground
				   const vector<float>& D, //input - distance map
				   //IMPORTANT! - B and B0 need to be initialized to zero by the calling function
				   vector<unsigned char>& B, //output - nodule segmentation 
				   vector<unsigned char>& B0, //output - nodule segmentation after the 1st phase
				   const vector<int>& vseeds, //input - seed positions in a set of triplets (x, y, z)
				   bool bConvexHull, //enable/disable a convex hull operation
				   bool bRefineSeg, //enable/disable an operation that refines segmentation boundary
				   int ndim,
				   const int* dims);


#endif /* ___SZ_EXTRACT_NODULE_H____ */