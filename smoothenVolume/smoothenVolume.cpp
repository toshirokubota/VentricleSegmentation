// smoothen.cpp : Defines the exported functions for the DLL application.
//

// Segmentor.cpp : Defines the entry point for the DLL application.
//

#include <iostream>
#include <algorithm>
using namespace std;
#include <stdio.h>

#ifdef MEX_DLL
#include <mex.h>
#include "mexFileIO.h"
#endif

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMiscOperations.h>

#include <szMyNeighborOp.h>
//#include <mytype.h>
#include <smoothen.h>
#include <SCpoint.h>

#ifdef MEX_DLL
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1)
	{
		mexErrMsgTxt("Usage: vargout = LCSmoothen(A, [niter tau lambda])");
		return;
	}

	//load image data
	int ndimA;
	const int* dimsA;
	mxClassID classA;
	vector<float> A;
	LoadData(A,prhs[0],classA,ndimA,&dimsA);

	int niter = 5;
	if(nrhs>=2) 
	{
		mxClassID classN;
		ReadScalar(niter, prhs[1], classN);
	}

	real tau = 1;
	if(nrhs>=3) 
	{
		mxClassID classN;
		ReadScalar(tau, prhs[2], classN);
	}
	
	real lambda = 1;
	if(nrhs>=4) 
	{
		mxClassID classN;
		ReadScalar(lambda, prhs[3], classN);
	}
	printf("LCSmoothen: niter = %d, tau = %f, lambda = %f\n", niter, tau, lambda);

#else
bool
liverCancerSegmentor(vector<float>& S,
					 const vector<unsigned short>& A,
					 real tau, real lambda, int niter,
					 int ndimA,
					 const int* dimsA) 
{
#endif
	int nvoxels = numberOfElements(ndimA,dimsA);

	vector<unsigned char> M(nvoxels, 1);
	//smoothen3D(A, M, niter, tau, lambda, dimsA);
	smoothen3L(A, M, niter, 0.5, dimsA);
	/*const int MatDimC = 10;
	float* gn = new float[MatDimC*MatDimC];
	float* gs = new float[MatDimC*MatDimC];
	float* gw = new float[MatDimC*MatDimC];
	float* ge = new float[MatDimC*MatDimC];
	if (ReadGammaMatricesC(gn, gw, ge, gs, tau, lambda, "GammaCube1-1b.txt") < 0)
	{
		mexErrMsgTxt("LCSegmentor: Failed to read gamma file.\n");
		return;
	}

	smoothen3(A, M, niter, gn, gs, gw, ge, tau, lambda, dimsA);*/

#ifdef MEX_DLL
	if(nlhs >= 1)
	{
		plhs[0] = StoreData(A,mxSINGLE_CLASS,ndimA,dimsA);
	}

	mexUnlock();
#else
	return bRet;
#endif
}


