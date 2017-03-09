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
#include <szFigureGroundSeparation.h>
#include <szExtractForeground.h>
#include <szMiscOperations.h>
#include <szFindSeeds.h>
#include <szConvexHull3D.h>
#include <szCutOutNodule.h>
#include <szParticle.h>
#include <szDistanceTransform.h>
#include <szConnectedComponent.h>
#include <szLocalExtrema.h>
#include <szTrace.h>
#include <szRetrace.h>
#include <szDistanceTransformNonIsotropic.h>

//Change these settings for another neighborhood system
const int NumNeighborsHere = NumNeighbors;
#define XOffsetHere XOffset
#define YOffsetHere YOffset
#define ZOffsetHere ZOffset

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

/*
onSurface3 with a general neighborhood
*/
template<class T>
bool
onSurface3B(const vector<T>& A, int x, int y, int z, const int* dims) 
{
	if(GetData3(A,x,y,z,dims[0],dims[1],dims[2],(T)0)) 
	{
		for(int m=0; m<NumNeighborsHere; ++m)
		{
			int x2 = x + XOffsetHere[m];
			int y2 = y + YOffsetHere[m];
			int z2 = z + ZOffsetHere[m];
			if(GetData3(A, x2, y2, z2, dims[0], dims[1], dims[2], (T)1)==0)
			{
				return true;
			}
		}
	}
	return false;
}

void
ExtactSurface(vector<unsigned char>& T, 
			  const vector<int>& C, 
			  const vector<unsigned char>& L, 
			  const int* dims)
{
	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<dims[1]; ++j)
		{
			for(int k=0; k<dims[0]; ++k)
			{
				if(onSurface3B(L, k, j, i, dims))
				{
					int cval = GetData3(C, k, j, i, dims[0], dims[1], dims[2], (int)0);
					if(cval)
					{
						SetData3(T, k, j, i, dims[0], dims[1], dims[2], (unsigned char)1);
					}
				}
			}
		}
	}
}

vector<CParticle>
collectSphericalFront(vector<unsigned char>& F,
					  const vector<unsigned char>& S, //foreground
					  const vector<float>& D,
					  const vector<CParticle>& vp, 
					  const vector<float>& v,
					  const int* dims)
{
	vector<unsigned char> L(S.size(), 0);
	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<dims[1]; ++j)
		{
			for(int k=0; k<dims[0]; ++k)
			{
				for(int n=0; n<vp.size(); ++n)
				{
					int x = vp[n].m_X;
					int y = vp[n].m_Y;
					int z = vp[n].m_Z;
					float dval = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
					float dx = (k-x)*v[0];
					float dy = (j-y)*v[1];
					float dz = (i-z)*v[2];
					if(dx*dx + dy*dy + dz*dz < dval*dval)
					{
						SetData3(L, k, j, i, dims[0], dims[1], dims[2], (unsigned char) 1);
					}
				}
			}
		}
	}
	vector<CParticle> vfront;
	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<dims[1]; ++j)
		{
			for(int k=0; k<dims[0]; ++k)
			{
				if(GetData3(L, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					int sval = GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0);
					if(sval)
					{
						SetData3(F, k, j, i, dims[0], dims[1], dims[2], (unsigned char) 1);
						if(onSurface3B(L, k, j, i, dims))
						{
							vfront.push_back(CParticle(k, j, i));
						}
					}
				}
			}
		}
	}
	return vfront;
}

/*
This version takes the voxel size into account.
*/
int
FreeGrowthHere(vector<int>& B,
		   const vector<float>& D,
		   const vector<CParticle>& vSeeds,
		   const vector<float>& v,
		   const int* dims)
{
	vector<CParticle> vfront = vSeeds;
	int cnt = 1;
	for(int i=0; i<vfront.size(); ++i)
	{
		vfront[i].m_Life = 0;
		SetVoxel(B, vfront[i], cnt, dims);
	}

	float time = 0;
	vector<int> steps(3, 0); 
	while(!vfront.empty())
	{
		cnt++;
		vector<CParticle> vfront2 = vfront;
		vector<CParticle> vnew;

		//compute the next time point
		float timeNext;
		for(int i=0; i<steps.size(); ++i)
		{
			float t = (steps[i]+1)*v[i];
			if(i==0 || t < timeNext)
			{
				timeNext = t;
			}
		}

		time = timeNext;
		for(int i=0; i<steps.size(); ++i)
		{
			float t = (steps[i]+1)*v[i];
			if(t <= time)
			{
				steps[i]++;
			}
		}	

		for(int i=0; i<vfront.size(); ++i)
		{
			int x1=vfront[i].m_X;
			int y1=vfront[i].m_Y;
			int z1=vfront[i].m_Z;
			float dval1 = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
			bool bDone = true;
			for(int j=0; j<NumNeighborsHere; ++j)
			{
				int x2 = x1+XOffsetHere[j];
				int y2 = y1+YOffsetHere[j];
				int z2 = z1+ZOffsetHere[j];
				if(!GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (int)1))
				{
					float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
					if(dval2>0)
					{
						if(dval1 > dval2)
						{
							float len2 = sqrt((float)v[0]*v[0]*XOffsetHere[j]*XOffsetHere[j] + v[1]*v[1]*YOffsetHere[j]*YOffsetHere[j] + v[2]*v[2]*ZOffsetHere[j]*ZOffsetHere[j]);
							if(time - vfront[i].m_Life >= len2)
							{
								CParticle pn(x2, y2, z2, time);
								if(find(vnew.begin(), vnew.end(), pn) == vnew.end() && find(vfront.begin(), vfront.end(), pn) == vfront.end())
								{
									vnew.push_back(pn);
								}
							}
							else
							{
								bDone = false;
							}
						}
					}
				}
			}
			if(bDone)
			{
				vfront2.erase(find(vfront2.begin(), vfront2.end(), vfront[i]));
			}
		}
		for(int i=0; i<vnew.size(); ++i)
		{
			SetVoxel(B, vnew[i], cnt, dims);
		}
		vfront2.insert(vfront2.end(), vnew.begin(), vnew.end());
		vfront = vfront2;
	}
	return cnt;
}


/*
Each surface voxel (T) in the free-trace volume (S) has to come from steepest descent.
Hence, the voxel in the neighborhood with the highest ascent has to have the trace sequence
number of one less than the surface voxel in question.
*/
void
CheckTheGrowth(vector<unsigned char>& U,
			   const vector<unsigned char>& T, 
			   const vector<int>& S, 
			   const vector<float>& D, 
			   const int* dims)
{
	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<dims[1]; ++j)
		{
			for(int k=0; k<dims[0]; ++k)
			{
				unsigned char tval = GetData3(T, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0);
				if(tval)
				{
					int sval = GetData3(S, k, j, i, dims[0], dims[1], dims[2], (int)0);
					if(sval == 1)
					{
						//this is within the initial sphere.  Accept it.
						SetData3(U, k, j, i, dims[0], dims[1], dims[2], tval);
						continue;
					}
					float dval = GetData3(D, k, j, i, dims[0], dims[1], dims[2], (float)0);
					float maxrate = 0;
					vector<CParticle> vsource;
					for(int m=0; m<NumNeighborsHere; ++m)
					{
						int x2 = k + XOffsetHere[m];
						int y2 = j + YOffsetHere[m];
						int z2 = i + ZOffsetHere[m];
						float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
						float len = sqrt((float)XOffsetHere[m]*XOffsetHere[m] + YOffsetHere[m]*YOffsetHere[m] + ZOffsetHere[m]*ZOffsetHere[m]);
						float rate = (dval2 - dval)/len;
						if(rate > maxrate)
						{
							maxrate = rate;
							vsource.clear();
							vsource.push_back(CParticle(x2, y2, z2));
						}
						else if(rate == maxrate)
						{
							vsource.push_back(CParticle(x2, y2, z2));
						}
					}
					bool bSuccess = false;
					for(int n=0; n<vsource.size(); ++n)
					{
						int sval2 = GetVoxel(S, vsource[n], (int)0, dims);
						if(sval2 == sval - 1)
						{
							bSuccess = true;
							break;
						}
					}
					if(bSuccess)
					{
						SetData3(U, k, j, i, dims[0], dims[1], dims[2], tval);
					}
				}
			}
		}
	}
}

/*
For each connected component of T, the region should grow gradually, registering all growth numbers in sequence.
When the growh enters non-convex part of the foreground, some growth will not meet the first criterion (CheckTheGrowth()).
Thus, the connected componet will be broken into multiple ones in U.  Those components in U that do not have the minimum
growth numbers are resulted from this break-up, and are considered spurious.  We thus delete them.
*/
void
CheckTheGrowth2(vector<unsigned char>& V,
				const vector<unsigned char>& U, 
				const vector<unsigned char>& T, 
				const vector<int>& S, 
				const int* dims)
{
	vector<int> C(S.size(), 0);
	vector<int> nbh = MakeEightNeighborhood(3,dims);
	int nc = ConnectedComponentAnalysisBigger(C, T, nbh, (unsigned char)0, 3, dims); 
	printf("There are %d clusters in T.\n", nc);

	for(int n=1; n<=nc; ++n)
	{
		vector<unsigned char> L(U.size(), 0);
		int minLabel = 0;
		for(int i=0; i<L.size(); ++i)
		{
			if(C[i]==n && U[i]) ///Including U[i] correct?
			{
				if(minLabel == 0)
				{
					minLabel = S[i];
				}
				else if(minLabel > S[i])
				{
					minLabel = S[i];
				}
				L[i]=1;
			}
		}
		vector<int> D(S.size(), 0);
		int nc2 = ConnectedComponentAnalysisBigger(D, L, nbh, (unsigned char)0, 3, dims); 
		printf("minlabel = %d.  There are %d clusters in Component %d.\n", minLabel, nc2, n);

		vector<bool> bFlag(nc2+1, false);
		for(int i=0; i<D.size(); ++i)
		{
			if(D[i]>0)
			{
				if(S[i] == minLabel)
				{
					bFlag[D[i]]=true;
				}
			}
		}
		for(int i=1; i<bFlag.size(); ++i)
		{
			if(bFlag[i])
			{
				printf("Component %d is accepted.\n", i);
			}
		}
		for(int i=0; i<D.size(); ++i)
		{
			if(D[i]>0 && bFlag[D[i]])
			{
				V[i] = 1;
			}
		}
	}
}

#ifdef MEX_DLL
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 3 || nlhs < 0)
	{
		//mexErrMsgTxt("Usage: [C [S D]] = TraceTest(S, L, seed)");
		mexErrMsgTxt("Usage: [C] = TraceTest(L, v)");
		return;
	}

	//load figure-ground segmentation
	int ndimL;
	const int* dimsL;
	mxClassID classL;
	vector<unsigned char> L;
	LoadData(L,prhs[0],classL,ndimL,&dimsL);

	int ndimv;
	const int* dimsv;
	mxClassID classv;
	vector<float> v;
	LoadData(v,prhs[1],classv,ndimv,&dimsv);

	vector<int> vseed;
	int ndimT;
	const int* dimsT;
	mxClassID classT;
	LoadData(vseed,prhs[2],classT,ndimT,&dimsT);

	int nvoxels = numberOfElements(ndimL,dimsL);
#else
bool
DmapTrace(vector<unsigned char>& C, //segmentation result
		  vector<unsigned char>& B, //free-trace result
		  const vector<unsigned char>& A,
		  const vector<int>& vseed,		  
		  int ndimA,
		  const int* dimsA
		  ) 
{
	int nvoxels = numberOfElements(ndimA,dimsA);
#endif

	vector<CParticle> vpclick; //(1, CParticle(vseed[0], vseed[1], vseed[2]));
	for(int i=0; i<vseed.size(); i+=3)
	{
		printf("straightGrowth: Seed = (%d, %d, %d)\n", vseed[i], vseed[i+1], vseed[i+2]);
		vpclick.push_back(CParticle(vseed[i], vseed[i+1], vseed[i+2]));
	}

	vector<float> D(nvoxels,0);
	vector<unsigned char> iL(nvoxels,0);
	for(int i=0; i<nvoxels; ++i)
	{
		iL[i] = L[i] ? 0: 1;
	}
	DistanceTransformEuclidF(D, iL, v, ndimL, dimsL);

	vector<unsigned char> F(nvoxels, 0);
	vector<CParticle> vfront = collectSphericalFront(F, L, D, vpclick, v, dimsL);

	printf("There are %d front voxels.\n", vfront.size());
	vector<int> S(nvoxels, 0);
	for(int i=0; i<nvoxels; ++i)
	{
		if(F[i]) S[i] = 1;
	}
	FreeGrowthHere(S, D, vfront, v, dimsL);		

	vector<unsigned char> T(nvoxels, 0);
	ExtactSurface(T, S, L, dimsL);

	vector<unsigned char> T2(T.size(), 0);
	CheckTheGrowth(T2, T, S, D, dimsL);

	vector<unsigned char> T3(T.size(), 0);
	CheckTheGrowth2(T3, T2, T, S, dimsL);

#ifdef MEX_DLL
	if(nlhs >= 1)
	{
		plhs[0] = StoreData(S,mxINT32_CLASS,ndimL,dimsL);
	}
	if(nlhs >= 2)
	{
		plhs[1] = StoreData(T3,mxUINT8_CLASS,ndimL,dimsL);
	}
	if(nlhs >= 3)
	{
		plhs[2] = StoreData(T2,mxUINT8_CLASS,ndimL,dimsL);
	}
	if(nlhs >= 4)
	{
		plhs[3] = StoreData(T,mxUINT8_CLASS,ndimL,dimsL);
	}
	if(nlhs >= 5)
	{
		plhs[4] = StoreData(D,mxSINGLE_CLASS,ndimL,dimsL);
	}
	mexUnlock();
#else
	return true;
#endif
}

