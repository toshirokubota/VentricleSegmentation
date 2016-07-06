#ifndef ___NODULE_SEGMENT_FIGURE_GROUND_SEPARATION_H___
#define ___NODULE_SEGMENT_FIGURE_GROUND_SEPARATION_H___

#include <vector>
using namespace std;

bool
figureGroundSeparation(vector<unsigned char>& L,
					   const vector<unsigned short>& A,
					   const vector<float>& v,
					   int noduleType, 
					   int ndimA,
					   const int* dimsA);

#endif /* ___NODULE_SEGMENT_FIGURE_GROUND_SEPARATION_H___ */
