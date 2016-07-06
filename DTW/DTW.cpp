#include <iostream>
using namespace std;
#include <stdio.h>

#ifdef MEX_DLL
#include <mex.h>
#include "mexFileIO.h"
#endif

#include <vector>
#include <algorithm>
using namespace std;
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szDefaultParam.h>
#include <szMiscOperations.h>
#include <szParticle.h>
#include <DynamicTimeWarp.h>


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 2 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [z cost] = DTW(x, y, params)");
		return;
	}

	//load an input signal 
	int ndimX;
	const int* dimsX;
	mxClassID classX;
	vector<float> x;
	LoadData(x,prhs[0],classX,ndimX,&dimsX);

	//load a reference signal (possibly multiple of them). 
	int ndimY;
	const int* dimsY;
	mxClassID classY;
	vector<float> Y;
	LoadData(Y, prhs[1], classX, ndimX, &dimsX);

	vector<float> params;
	int ndimP;
	const int* dimsP;
	mxClassID classP;
	LoadData(params,prhs[2],classP,ndimP,&dimsP);

	float minCost = DynamicTimeWarp(x, Y, params);

	if(nlhs >= 1)
	{
		plhs[0] = StoreScalar(minCost,mxSINGLE_CLASS);
	}
	mexUnlock();
}

