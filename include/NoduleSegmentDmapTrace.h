#ifndef ___NODULE_SEGMENT_DMAP_TRACE_H___
#define ___NODULE_SEGMENT_DMAP_TRACE_H___

#include <vector>
using namespace std;

bool
DmapTrace(vector<unsigned char>& C, //segmentation result
		  vector<unsigned char>& B, //free-trace result
		  const vector<unsigned char>& A, //foreground
		  const vector<int>& vseed,		  
		  int ndimA,
		  const int* dimsA
		  );

#endif /* ___NODULE_SEGMENT_DMAP_TRACE_H___ */
