#ifndef ___SZ_TRACE_H___
#define ___SZ_TRACE_H___
#include <vector>
using namespace std;
#include <szParticle.h>

void
FreeGrowth(vector<unsigned char>& B,
		   const vector<float>& D,
		   const vector<CParticle>& vSeeds,
		   const int* dims);


void
ConstrainedGrowth(vector<unsigned char>& B,
				  const vector<float>& D,
				  vector<unsigned char>& P,
				  int x, int y, int z,
				  const int* dims);

#endif /* ___SZ_TRACE_H___ */
