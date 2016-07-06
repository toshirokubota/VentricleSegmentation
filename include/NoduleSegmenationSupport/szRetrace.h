#ifndef ___SZ_RETRACE_H___
#define ___SZ_RETRACE_H___
#include <vector>
using namespace std;
#include <szParticle.h>

int
traceUp(vector<int>& M, 
		const vector<float>& D, 
		const CParticle& p, 
		const vector<CParticle>& vppeaks, 
		const int* dims);

vector<CParticle>
collectTraced(const vector<int>& M,
			  const vector<CParticle>& vppeaks,
			  const int* dims);

#endif /* ___SZ_RETRACE_H___ */
