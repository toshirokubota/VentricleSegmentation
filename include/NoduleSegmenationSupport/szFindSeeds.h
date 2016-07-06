#ifndef _SZ_FIND_SEEDS_H_
#define _SZ_FIND_SEEDS_H_

#include <vector>
using namespace std;
#include <szParticle.h>

vector<int>
FindSeedsThres(const vector<unsigned char>& C,
          vector<float>& D,
          vector<float>& Sp,
          vector<unsigned char>& P,
		  vector<float>& vSpherical,
          int w,
		  int wcore, 
          float selection_threshold,
          const int* dims);

vector<int>
FindSeedsThresNew(const vector<unsigned char>& C,
          vector<float>& D,
          vector<float>& Sp,
          vector<unsigned char>& P,
		  vector<float>& vSpherical,
          int w,
		  int wcore, 
          float selection_threshold,
		  int numseeds,
          const int* dims);

bool traceUphill(const vector<float>& D,
				 int x0, int y0, int z0,
				 int x, int y, int z,
				 int& xf, int& yf, int& zf,
				 float dthres,
				 const int* dims);

bool traceUphillAll(const vector<float>& D,
					vector<unsigned char>& M,
					int x, int y, int z,
					vector<CParticle>& vppeak,
					float dthres,
					const int* dims);

vector<int>
refineSeeds(const vector<int>& vseeds,
			const vector<float>& D,
			const vector<unsigned char>& L,
			float dval, 
			const int* dims);

vector<CParticle>
FindSeedsThresNew(const vector<unsigned char>& L,  //INPUT - foreground segmentation
          vector<float>& D,          //OUTPUT - distance map
          vector<float>& Sp,          //OUTPUT - sphericity map
          vector<unsigned char>& P,  //OUTPUT - local maxima of the sphericity map
		  vector<float>& vSpherical, //OUTPUT - spherical values of seeds
		  vector<CParticle>& vclick, //INPUT - click points
          int w,
          float selection_threshold,
          const int* dims);

vector<CParticle>
FindSeedsThresV3(const vector<float>& Sp,          //INPUT - sphericity map
				 const vector<unsigned char>& L,          //INPUT - foreground map
				 const vector<float>& D,          //INPUT - distance map
				 const vector<CParticle>& vclick, //INPUT - click points
				 const vector<float>& v,
				 float search_core,
				 float perc,
				 const int* dims);

#endif /* _SZ_FIND_SEEDS_H_ */