#ifndef ___NODULE_SEGMENT_FIND_SEEDS_H___
#define ___NODULE_SEGMENT_FIND_SEEDS_H___

#include <vector>
using namespace std;

vector<int>
findSeeds(const vector<unsigned short>& A, //original data
		  const vector<unsigned char>& L,  //foreground
		  vector<float>& D,	//distance map of the foreground
		  vector<float>& Sp,//spericity map
		  vector<unsigned char> P, //local maximum of sphericity values
		  int findSeedWidth,
		  int findSeedCoreWidth,
		  float findSeedThres,
		  int ndimA,
		  const int* dimsA);

#endif /* ___NODULE_SEGMENT_FIND_SEEDS_H___ */
