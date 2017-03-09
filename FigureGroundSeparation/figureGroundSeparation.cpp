#include <iostream>
using namespace std;
#include <stdio.h>

#ifdef MEX_DLL
#include <mex.h>
#include "mexFileIO.h"
#endif

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szReactionDiffusion.h>
#include <szMiscOperations.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 2 || nlhs < 1)
	{
		mexErrMsgTxt("Usage: vargout = FigureGroundSeparation(A, v, [modes, numIter, lambda])");
		return;
	}

	//load image data
	int ndimA;
	const int* dimsA;
	mxClassID classA;
	vector<float> Q;
	LoadData(Q,prhs[0],classA,ndimA,&dimsA);

	//load voxel size
	int ndimv;
	const int* dimsv;
	mxClassID classv;
	vector<float> v;
	LoadData(v,prhs[1],classv,ndimv,&dimsv);
	if(v.size() < 3)
	{
		mexErrMsgTxt("FigureGroundSeparation: The voxel size needs three elements.");
		return;
	}
	printf("voxel size: %f, %f, %f\n", v[0], v[1], v[2]);

	int numiter = DefaultRDNumIteration;
	if(nrhs>=3) 
	{
		mxClassID classMode;
		ReadScalar(numiter,prhs[2],classMode);
	}
	//float lambda = DefaultRDLambda;
	vector<float> coeff;
	if(nrhs>=4) 
	{
		mxClassID classC;
		int ndimC;
		const int* dimsC;
		LoadData(coeff,prhs[3],classC,ndimC,&dimsC);
	}
	else
	{
		coeff.push_back(1); coeff.push_back(0); coeff.push_back(0); coeff.push_back(1); 
	}
	float dfrate = 0.125;
	if(nrhs>=5) 
	{
		mxClassID classMode;
		ReadScalar(dfrate,prhs[4],classMode);
	}
	float rate = 0.25;
	if(nrhs>=6) 
	{
		mxClassID classMode;
		ReadScalar(rate,prhs[5],classMode);
	}
	int minSize = 300;
	if (nrhs >= 7)
	{
		mxClassID classMode;
		ReadScalar(minSize, prhs[6], classMode);
	}

	int nvoxels = numberOfElements(ndimA, dimsA);
	vector<unsigned char> L(nvoxels,(float)0);

	printf("numIter=%d, coeff[0]=%f, coeff[3]=%f, rate=%f, fillsize=%d\n", numiter, coeff[0], coeff[3], rate, minSize);
	if (ndimA == 3)
	{
		float dfspeed[] = { dfrate / v[2], dfrate / v[2], dfrate / v[2], dfrate / v[2], dfrate / v[2], dfrate / v[2] }; //{1,1,1,1,0.5,0.5};
		//printf("dfspeed: %f, %f, %f, %f, %f, %f\n", dfspeed[0], dfspeed[1], dfspeed[2], dfspeed[3], dfspeed[4], dfspeed[5]);

		for (int iter = 0; iter < numiter; ++iter) {
			ReactionDiffusion3D(Q, dimsA[0], dimsA[1], dimsA[2], coeff, dfspeed, rate);
		}
	}
	else if (ndimA == 4)
	{
		float dfspeed[] = { dfrate, dfrate, dfrate, dfrate, dfrate, dfrate, dfrate,  dfrate }; //{1,1,1,1,0.5,0.5};
		for (int iter = 0; iter < numiter; ++iter) {
			ReactionDiffusion4D(Q, dimsA[0], dimsA[1], dimsA[2], dimsA[3], coeff, dfspeed, rate);
		}
	}
	for(int i=0; i<Q.size(); ++i) {
		if(Q[i]>.5) L[i]=ForegroundColor;
		else L[i]=0;
	}

	FillBackgroundHoles(L, minSize, ForegroundColor, ndimA, dimsA);

#ifdef MEX_DLL
	if(nlhs >= 1)
	{
		plhs[0] = StoreData(L,mxUINT8_CLASS,ndimA,dimsA);
	}
	if(nlhs >= 2)
	{
		plhs[1] = StoreData(Q,mxSINGLE_CLASS,ndimA,dimsA);
	}	
	mexUnlock();
#else
	return bRet;
#endif
}

