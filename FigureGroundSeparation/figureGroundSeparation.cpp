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
#include <szFigureGroundSeparation.h>
#include <szMiscOperations.h>

#ifdef MEX_DLL
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
	vector<float> A;
	LoadData(A,prhs[0],classA,ndimA,&dimsA);

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

	int noduleType;
	if(nrhs>=3) 
	{
		mxClassID classT;
		ReadScalar(noduleType,prhs[2],classT);
	}
	else
	{
		noduleType = SZ_LABEL_NON_BACKGROUND;
	}
	int numiter = DefaultRDNumIteration;
	if(nrhs>=4) 
	{
		mxClassID classMode;
		ReadScalar(numiter,prhs[3],classMode);
	}
	//float lambda = DefaultRDLambda;
	vector<float> coeff;
	if(nrhs>=5) 
	{
		mxClassID classC;
		int ndimC;
		const int* dimsC;
		LoadData(coeff,prhs[4],classC,ndimC,&dimsC);
	}
	else
	{
		coeff.push_back(1); coeff.push_back(0); coeff.push_back(0); coeff.push_back(1); 
	}
	float dfrate = 0.125;
	if(nrhs>=6) 
	{
		mxClassID classMode;
		ReadScalar(dfrate,prhs[5],classMode);
	}
	float rate = 0.25;
	if(nrhs>=7) 
	{
		mxClassID classMode;
		ReadScalar(rate,prhs[6],classMode);
	}
	int nvoxels = numberOfElements(ndimA,dimsA);
	vector<unsigned char> L(nvoxels,(float)0);

	printf("mode = %d, numIter=%d, coeff[0]=%f, coeff[3]=%f, rate=%f\n", noduleType, numiter, coeff[0], coeff[3], rate);

#else
bool
figureGroundSeparation(vector<unsigned char>& L,
					   const vector<float>& A,
					   const vector<float>& v,
					   int noduleType, 
					   int ndimA,
					   const int* dimsA) 
{
	//int noduleType = SZ_LABEL_NON_BACKGROUND;
	int numiter = DefaultRDNumIteration;
	float lambda = DefaultRDLambda;
	int nvoxels = numberOfElements(ndimA,dimsA);
#endif
	int i;
	vector<float> Q(nvoxels,(float)0);

	/*bool bRet = FigureGroundSeparation(
		A,
		L,
		Q,
		noduleType,
		numiter,
		lambda,
		dimsA);*/

	float maxA = A[0], minA = A[0];
	for(int i=1; i<A.size(); ++i)
	{
		if(A[i] > maxA) maxA = A[i]; 
		if(A[i] < minA) minA = A[i];
	}

	if(maxA <= 1.0 && minA >= 0) //the input is a probability volume
	{
		printf("FigureGroundSeparation: the input volume is interpreted as a probabilistic one.\n");
		for(int i=0; i<A.size(); ++i)
		{
			Q[i] = A[i];
		}
	}
	else
	{
		printf("FigureGroundSeparation: the input volume is interpreted as a raw-valued one.\n");
		ReactionDiffusionInitialize(A, Q, noduleType, dimsA);
	}

	float minV = dfrate*min(v[0], min(v[1], v[2]));
	//float dfspeed[]={1.0/v[2], 1.0/v[2], 1.0/v[2], 1.0/v[2], 1.0/v[2], 1.0/v[2]}; 
	//float dfspeed[]= {minV/v[0],minV/v[0],minV/v[1], minV/v[1], minV/v[2], minV/v[2]}; //{1,1,1,1,0.5,0.5};
	float dfspeed[]= {dfrate/v[2], dfrate/v[2], dfrate/v[2], dfrate/v[2], dfrate/v[2], dfrate/v[2]}; //{1,1,1,1,0.5,0.5};
	//float dfspeed[]= {dfrate/v[0], dfrate/v[0], dfrate/v[1], dfrate/v[1], dfrate/v[2], dfrate/v[2]}; //{1,1,1,1,0.5,0.5};
	//float dfspeed[]= {1.0/16.0, 1.0/16.0, 1.0/16.0, 1.0/16.0, .125, .125};
	printf("dfspeed: %f, %f, %f, %f, %f, %f\n", dfspeed[0], dfspeed[1], dfspeed[2], dfspeed[3], dfspeed[4], dfspeed[5]);
	float reacWeight=1.0;

	for(int iter=0; iter<numiter; ++iter) {
		ReactionDiffusionNew(Q,dimsA[0],dimsA[1],dimsA[2],coeff, dfspeed, rate);
	}

	for(i=0; i<Q.size(); ++i) {
		if(Q[i]>.5) L[i]=ForegroundColor;
		else L[i]=0;
	}

	int minSize = 300;
	FillBackgroundHoles(L, minSize, ForegroundColor, dimsA);

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

