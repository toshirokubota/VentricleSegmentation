#ifndef ___SZ_CUT_OUT_NODULE_H___
#define ___SZ_CUT_OUT_NODULE_H___

#include <vector>
using namespace std;

#include "szParticle.h"

void
preConvexHull(vector<unsigned char>& S,   //OUTPUT - must be initialized to zero
              const vector<unsigned char>& C,   //INPUT - segmentation map
              const int* dims);

void
postConvexHull(vector<unsigned char>& B,    //OUTPUT - must be initialized to zero
               const vector<unsigned char>& C,    //INPUT - convex-hull map
               const vector<unsigned char>& L,   //INPUT - foreground segmentation map
               const int* dims);             //volume size

vector<CParticle>
TraceNew(vector<unsigned char>& B,
		 const vector<unsigned char>& L, 
		 const vector<float>& D, 
		 vector<CParticle> vseeds,
		 float scale,
		 bool (*decision_func)(float, float),
		 const int* dims);

vector<CParticle>
TraceNew2(vector<unsigned char>& B,
		 const vector<unsigned char>& L, 
		 const vector<float>& D, 
		 vector<CParticle> vseeds,
		 float scale,
		 bool (*decision_func)(float, float),
		 const int* dims);

void
MergeLabels(vector<unsigned char>& A,
            vector<unsigned char>& B,
            const int* dims);

bool
AlwaysTrace(float a, float b);

bool
StrictDownHillTrace(float a, float b);

bool
NonStrictDownHillTrace(float a, float b);

#endif /* ___SZ_CUT_OUT_NODULE_H___ */