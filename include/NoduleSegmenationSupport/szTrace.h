#ifndef ___SZ_TRACE_H___
#define ___SZ_TRACE_H___
#include <vector>
using namespace std;
#include <szParticle.h>

void
SampledTrace(vector<unsigned char>& B,
			 const vector<float>& D,
			 const vector<CParticle>& vSeeds,
			 const int* dims);

void
incrementalTrace(vector<unsigned char>& B,
				 const vector<float>& D,
				 const vector<CParticle>& vSeeds,
				 const int* dims);

int
FreeGrowth(vector<int>& B,
		   const vector<float>& D,
		   const vector<CParticle>& vSeeds,
		   bool bStrict,
		   const int* dims);

//6-neighborhood version of FreeGrowth
int
FreeGrowth6(vector<int>& B,
		   const vector<float>& D,
		   const vector<CParticle>& vSeeds,
		   bool bStrict,
		   const int* dims);

void
FreeGrowthAngle(vector<unsigned char>& B,
				const vector<float>& D,
				const vector<CParticle>& vseed,
				const int* dims);

void
ConstrainedGrowth(vector<unsigned char>& B,
				  const vector<float>& D,
				  vector<unsigned char>& P,
				  int x, int y, int z,
				  const int* dims);

void
NewRegionPartition(vector<unsigned char>& B,
				   const vector<float>& D,
				   vector<unsigned char>& P,
				   int x, int y, int z,
				   const int* dims);

void
YetNewRegionPartition(vector<unsigned char>& C,
					  const vector<unsigned char>& B, //free-growth
					  const vector<float>& D,	//distance map computed on the free-growh
					  vector<unsigned char>& P,	//local maxima of D
					  int x, int y, int z,		//seed point
					  const int* dims);

void
NegativeMigrateGrowth(vector<unsigned char>& B,
					  const vector<unsigned char>& S,
					  const vector<float>& D,
					  const vector<unsigned char>& P,
					  const vector<CParticle>& vseeds,
					  const int* dims);

#endif /* ___SZ_TRACE_H___ */
