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
		mexErrMsgTxt("Usage: [C] = DTWvolume(A, [params])");
		return;
	}

	//load an input volume 
	int ndim;
	const int* dims;
	mxClassID classX;
	vector<float> A;
	LoadData(A,prhs[0],classX,ndim,&dims);

	//load a reference signal (possibly multiple of them). 
	int ndimY;
	const int* dimsY;
	mxClassID classY;
	vector<float> y;
	LoadData(y, prhs[1], classX, ndimY, &dimsY);
	if (y.size() != dims[2])
	{
		mexErrMsgTxt("The length of a reference has to be the same with the 3rd dimension of the volume.");
		return;
	}

	vector<float> params;
	if (nrhs > 2)
	{
		int ndimP;
		const int* dimsP;
		mxClassID classP;
		LoadData(params, prhs[2], classP, ndimP, &dimsP);
	}

	vector<float> C(dims[0]*dims[1], 0.0f);
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			vector<float> x(dims[2]);
			for (int k = 0; k < dims[2]; ++k)
			{
				x[k] = GetData3(A, j, i, k, dims[0], dims[1], dims[2], 0.0f);
			}
			float minCost = DynamicTimeWarp(x, y, params);
			SetData2(C, j, i, dims[0], dims[1], minCost);
		}
	}
	if(nlhs >= 1)
	{
		const int dims2[] = { dims[0], dims[1] };
		plhs[0] = StoreData(C,mxSINGLE_CLASS, 2, dims2);
	}
	mexUnlock();
}

