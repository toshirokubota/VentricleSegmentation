#ifndef _SZ_FIND_SEEDS_H_
#define _SZ_FIND_SEEDS_H_

#include <vector>
using namespace std;

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

vector<int>
refineSeeds(const vector<int>& vseeds,
			const vector<float>& D,
			const vector<unsigned char>& L,
			float dval, 
			const int* dims);

vector<int>
AscendToBrightest(const vector<unsigned short>& A, 
				  const vector<float>& D, 
				  const vector<int>& vseeds, 
				  const int* dimsA);

#endif /* _SZ_FIND_SEEDS_H_ */