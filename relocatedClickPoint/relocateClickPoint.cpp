#include <iostream>
using namespace std;
#include <stdio.h>

#ifdef MEX_DLL
#include <mex.h>
#include "mexFileIO.h"
#endif

#include <vector>
#include <algorithm>
#include <set>
using namespace std;
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szDefaultParam.h>
#include <szMiscOperations.h>
#include <szParticle.h>
#include <szDistanceTransform.h>
#include <szDistanceTransformNonIsotropic.h>


template<class T>
bool
SetVoxel(vector<T>& A,
const CParticle& p,
T value,
const int* dims)
{
	return SetData3(A, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], value);
}

template<class T>
T
GetVoxel(const vector<T>& A,
const CParticle& p,
T defaultValue,
const int* dims)
{
	return GetData3(A, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], defaultValue);
}

vector<CParticle>
relocateClickPointsUpward(vector<float>& D, vector<CParticle>& vpclick, vector<float>& v, const int* dims)
{
	vector<CParticle> vpnew;
	for (int c = 0; c < vpclick.size(); ++c)
	{
		//vector<float> gr = Gradient(D, vpclick[c], dims);
		float dt = 1.0f;
		CParticle p = vpclick[c];
		float dval0 = GetVoxel(D, vpclick[c], 0.0f, dims);
		float dmax = dval0;
		vector<CParticle> Q(1, vpclick[c]);
		vector<CParticle> candiates(1, vpclick[c]);
		while (Q.empty() == false) {
			set<CParticle> S;
			float dmax0 = dmax;
			for (int i = 0; i < Q.size(); ++i)
			{
				CParticle p = Q[i];
				for (int iz = p.m_Z - 1; iz <= p.m_Z + 1; ++iz)
				{
					for (int iy = p.m_Y - 1; iy <= p.m_Y + 1; ++iy)
					{
						for (int ix = p.m_X - 1; ix <= p.m_X + 1; ++ix)
						{
							float dval = GetData3(D, ix, iy, iz, dims[0], dims[1], dims[2], 0.0f);
							if (dval > dmax)
							{
								dmax = dval;
								S.clear();
								S.insert(CParticle(ix, iy, iz));
							}
							else if (dval == dmax && dval > dmax0) 
							{
								S.insert(CParticle(ix, iy, iz));
							}
						}
					}
				}
			}
			if (S.empty())
			{
				vpnew.insert(vpnew.end(), Q.begin(), Q.end());
			}
			Q.clear();
			Q.insert(Q.end(), S.begin(), S.end());
		}
	}
	/*printf("Relocation results:\n");
	for (int i = 0; i < vpnew.size(); ++i)
	{
		printf("%d: (%d, %d, %d) with %f\n", i, vpnew[i].m_X, vpnew[i].m_Y, vpnew[i].m_Z, GetVoxel(D, vpnew[i], 0.0f, dims));
	}*/
	return vpnew;
}

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

	//load a click point. 
	int ndimY;
	const int* dimsY;
	mxClassID classY;
	vector<int> vp;
	LoadData(vp, prhs[1], classX, ndimY, &dimsY);
	vector<CParticle> vpclick; //(1, CParticle(vseed[0], vseed[1], vseed[2]));
	for (int i = 0; i<vp.size(); i += 3)
	{
		printf("relocateClickPoint: Seed = (%d, %d, %d)\n", vp[i], vp[i + 1], vp[i + 2]);
		vpclick.push_back(CParticle(vp[i], vp[i + 1], vp[i + 2]));
	}

	vector<float> voxelSize(3, 1.0f);

	int nvoxels = A.size();
	vector<float> D(nvoxels, 0);
	vector<unsigned char> iL(nvoxels, 0);
	for (int i = 0; i<nvoxels; ++i)
	{
		iL[i] = A[i] ? 0 : 1;
	}

	vector<float> vs(3, 1.0f);
	DistanceTransformEuclidF(D, iL, vs, ndim, dims);
	for (int i = 0; i<vpclick.size(); i++)
	{
		printf("Before: (%d, %d, %d),  %f\n", vpclick[i].m_X, vpclick[i].m_Y, vpclick[i].m_Z, GetVoxel(D, vpclick[i], -1.0f, dims));
	}
	vpclick = relocateClickPointsUpward(D, vpclick, voxelSize, dims);
	for (int i = 0; i<vpclick.size(); i++)
	{
		printf("After: (%d, %d, %d),  %f\n", vpclick[i].m_X, vpclick[i].m_Y, vpclick[i].m_Z, GetVoxel(D, vpclick[i], -1.0f, dims));
	}

	if(nlhs >= 1)
	{
		const int dims2[] = { vpclick.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < vpclick.size(); ++i)
		{
			SetData2(F, i, 0, dims2[0], dims2[1], vpclick[i].m_X);
			SetData2(F, i, 1, dims2[0], dims2[1], vpclick[i].m_Y);
			SetData2(F, i, 2, dims2[0], dims2[1], vpclick[i].m_Z);
		}
		plhs[0] = StoreData(F,mxINT32_CLASS, 2, dims2);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreData(D, mxSINGLE_CLASS, ndim, dims);
	}
	mexUnlock();
}

