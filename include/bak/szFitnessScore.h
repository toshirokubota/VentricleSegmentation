#ifndef _SZ_FITNESS_SCORE_H_
#define _SZ_FITNESS_SCORE_H_

#include <vector>
using namespace std;

int NumberOfFeatures();

vector<double>
fitnessScore(vector<unsigned char>& B, //non-convex region - before convex Hull
			 vector<float>& D,         //distance map
			 vector<float>& Sp,        //special sphericity measure
			 vector<float>& v,         //voxel size
			 int ndim,
			 const int* dims);

vector<double>
fitnessScoreNew(const vector<unsigned short>& A, //image
				const vector<unsigned char>& B, //nodule segmentation
				const vector<unsigned char>& L, //foreground segmentation
				const vector<float>& D,         //distance map
				const vector<float>& Sp,         //sphericity
				const vector<unsigned char>& P,  //sphericity local maxima
				const vector<int>& vseeds,      //seed positions
				const vector<float>& v,         //voxel size
				int ndim,
				const int* dims);

#endif /* _SZ_FITNESS_SCORE_H_ */