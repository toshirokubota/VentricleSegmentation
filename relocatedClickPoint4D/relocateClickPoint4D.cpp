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
#include <szParticle4D.h>
#include <szDistanceTransform.h>
#include <szDistanceTransformNonIsotropic.h>

template<class T>
bool
SetVoxel(vector<T>& A,
const CParticle4D& p,
T value,
const int* dims)
{
	return SetData4(A, p.m_X, p.m_Y, p.m_Z, p.m_T, dims[0], dims[1], dims[2], dims[3], value);
}

template<class T>
T
GetVoxel(const vector<T>& A,
const CParticle4D& p,
T defaultValue,
const int* dims)
{
	return GetData4(A, p.m_X, p.m_Y, p.m_Z, p.m_T, dims[0], dims[1], dims[2], dims[3], defaultValue);
}

vector<CParticle4D>
relocateClickPointsUpward(vector<float>& D, vector<CParticle4D>& vpclick, vector<float>& v, const int* dims)
{
	vector<CParticle4D> vpnew;
	for (int c = 0; c < vpclick.size(); ++c)
	{
		//vector<float> gr = Gradient(D, vpclick[c], dims);
		float dt = 1.0f;
		//CParticle p = vpclick[c];
		float dval0 = GetVoxel(D, vpclick[c], 0.0f, dims);
		float dmax = dval0;
		vector<CParticle4D> Q(1, vpclick[c]);
		vector<CParticle4D> candiates(1, vpclick[c]);
		while (Q.empty() == false) {
			set<CParticle4D> S;
			float dmax0 = dmax;
			for (int i = 0; i < Q.size(); ++i)
			{
				CParticle4D p = Q[i];
				for (int it = p.m_T - 1; it <= p.m_T + 1; ++it)
				{
					for (int iz = p.m_Z - 1; iz <= p.m_Z + 1; ++iz)
					{
						for (int iy = p.m_Y - 2; iy <= p.m_Y + 2; ++iy)
						{
							for (int ix = p.m_X - 2; ix <= p.m_X + 2; ++ix)
							{
								float dval = GetData4(D, ix, iy, iz, it, dims[0], dims[1], dims[2], dims[3], 0.0f);
								if (dval > dmax)
								{
									dmax = dval;
									S.clear();
									S.insert(CParticle4D(ix, iy, iz, it));
								}
								else if (dval == dmax && dval > dmax0)
								{
									S.insert(CParticle4D(ix, iy, iz, it));
								}
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
	vector<CParticle4D> vpclick; //(1, CParticle(vseed[0], vseed[1], vseed[2]));
	for (int i = 0; i<vp.size(); i += 4)
	{
		printf("relocateClickPoint: Seed = (%d, %d, %d)\n", vp[i], vp[i + 1], vp[i + 2]);
		vpclick.push_back(CParticle4D(vp[i], vp[i + 1], vp[i + 2],  vp[i+3]));
	}

	vector<float> voxelSize(4, 1.0f);

	int nvoxels = A.size();
	vector<float> D(nvoxels, 0);
	vector<unsigned char> iL(nvoxels, 0);
	for (int i = 0; i<nvoxels; ++i)
	{
		iL[i] = A[i] ? 0 : 1;
	}

	vector<float> vs(4, 1.0f);
	DistanceTransformEuclidF(D, iL, vs, ndim, dims);
	for (int i = 0; i<vpclick.size(); i++)
	{
		printf("Before: (%d, %d, %d, %d),  %f\n", 
			vpclick[i].m_X, vpclick[i].m_Y, vpclick[i].m_Z, vpclick[i].m_T, GetVoxel(D, vpclick[i], -1.0f, dims));
	}
	vpclick = relocateClickPointsUpward(D, vpclick, voxelSize, dims);
	for (int i = 0; i<vpclick.size(); i++)
	{
		printf("After: (%d, %d, %d, %d),  %f\n", 
			vpclick[i].m_X, vpclick[i].m_Y, vpclick[i].m_Z, vpclick[i].m_T, GetVoxel(D, vpclick[i], -1.0f, dims));
	}

	if(nlhs >= 1)
	{
		const int dims2[] = { vpclick.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < vpclick.size(); ++i)
		{
			SetData2(F, i, 0, dims2[0], dims2[1], vpclick[i].m_X);
			SetData2(F, i, 1, dims2[0], dims2[1], vpclick[i].m_Y);
			SetData2(F, i, 2, dims2[0], dims2[1], vpclick[i].m_Z);
			SetData2(F, i, 3, dims2[0], dims2[1], vpclick[i].m_T);
		}
		plhs[0] = StoreData(F,mxINT32_CLASS, 2, dims2);
	}
	if (nlhs >= 2)
	{
		plhs[1] = StoreData(D, mxSINGLE_CLASS, ndim, dims);
	}
	mexUnlock();
}

