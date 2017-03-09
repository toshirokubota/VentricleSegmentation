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
#include <map>
using namespace std;
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szDefaultParam.h>
//#include <szFigureGroundSeparation.h>
#include <szExtractForeground.h>
#include <szMiscOperations.h>
#include <szFindSeeds.h>
#include <szConvexHull3D.h>
#include <szCutOutNodule.h>
#include <szParticle4D.h>
#include <szDistanceTransform.h>
#include <szConnectedComponent.h>
#include <szLocalExtrema.h>
#include <szTrace.h>
#include <szRetrace.h>
#include <szDistanceTransformNonIsotropic.h>
//#define FOUR_NEIGHBORHOOD

#ifdef FOUR_NEIGHBORHOOD
/*int XOffsetArray[] = { -1, 1, 0, 0, 0, 0, 0, 0 };
int YOffsetArray[] = { 0, 0, -1, 1, 0, 0, 0, 0 };
int ZOffsetArray[] = { 0, 0, 0, 0, -1, 1, 0, 0 };
int TOffsetArray[] = { 0, 0, 0, 0, 0, 0, -1, 1 };*/
/*int XOffsetArray[] = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, 0, 0 };
int YOffsetArray[] = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1, 0, 0 };
int ZOffsetArray[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0 };
int TOffsetArray[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1 };*/
int XOffsetArray[] = { -1, 0, 1, -1, 1, -1, 0, 1, 0, 0, 0, 0 };
int YOffsetArray[] = { -1, -1, -1, 0, 0, 1, 1, 1, 0, 0, 0, 0 };
int ZOffsetArray[] = { 0, 0, 0, 0, 0, 0, 0, 0, -1, 1, 0, 0 };
int TOffsetArray[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 1 };
const int NumNeighborsHere = sizeof(XOffsetArray) / sizeof(XOffsetArray[0]);
#else
//Change these settings for another neighborhood system
const int NumNeighborsHere = 3 * NumNeighbors + 2;
int TOffsetArray[NumNeighborsHere];
int XOffsetArray[NumNeighborsHere];
int YOffsetArray[NumNeighborsHere];
int ZOffsetArray[NumNeighborsHere];
#endif

#define XOffsetHere XOffsetArray
#define YOffsetHere YOffsetArray
#define ZOffsetHere ZOffsetArray
#define TOffsetHere TOffsetArray

void
initializeOffsetArrays()
{
#ifdef FOUR_NEIGHBORHOOD
#else
	int idx = 0;
	for (int m = -1; m <= 1; ++m)
	{
		for (int i = -1; i <= 1; ++i)
		{
			for (int j = -1; j <= 1; ++j)
			{
				for (int k = -1; k <= 1; ++k)
				{
					if (m == 0 && i == 0 && j == 0 && k == 0)
					{
						continue;
					}
					else
					{
						XOffsetArray[idx] = k;
						YOffsetArray[idx] = j;
						ZOffsetArray[idx] = i;
						TOffsetArray[idx] = m;
						idx++;
					}
				}
			}
		}
	}
	for (int m = 0; m < NumNeighborsHere; ++m)
	{
		//printf("%d: %d %d %d %d\n", m, XOffsetHere[m], YOffsetHere[m], ZOffsetHere[m], TOffsetHere[m]);
	}
#endif
}

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

/*
onSurface3 with a general neighborhood
*/
template<class T>
bool
onSurface4(const vector<T>& A, int x, int y, int z, int t, const int* dims) 
{
	if(GetData4(A,x,y,z,t, dims[0],dims[1],dims[2],dims[3], (T)0)) 
	{
		for(int m=0; m<NumNeighborsHere; ++m)
		{
			int x2 = x + XOffsetHere[m];
			int y2 = y + YOffsetHere[m];
			int z2 = z + ZOffsetHere[m];
			int t2 = t + TOffsetHere[m];
			if (GetData4(A, x2, y2, z2, t2, dims[0], dims[1], dims[2], dims[3], (T)1) == 0)
			{
				return true;
			}
		}
	}
	return false;
}

void removeIsolated(vector<unsigned char>& F, 
	int minsize,
	const int* dims)
{
	vector<int> nbh = MakeEightNeighborhood(2, dims); //MakeEightNeighborhood(4, dims);
	int nremoved = 0;
	for (int m = 0; m<dims[3]; ++m)
	{
		for (int i = 0; i < dims[2]; ++i)
		{
			vector<unsigned char> L(dims[1] * dims[0]);
			for (int j = 0; j < dims[1]; ++j)
			{
				for (int k = 0; k < dims[0]; ++k)
				{
					unsigned char val = GetData4(F, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0);
					if (val)
					{
						SetData2(L, k, j, dims[0], dims[1], (unsigned char)1);
					}
				}
			}
			vector<int> C(L.size(), 0);
			int nc = ConnectedComponentAnalysisBigger(C, L, nbh, (unsigned char)0, 2, dims);
			vector<int> count(nc + 1, 0);
			for (int k = 0; k < C.size(); ++k)
			{
				count[C[k]]++;
			}

			for (int j = 0; j < dims[1]; ++j)
			{
				for (int k = 0; k < dims[0]; ++k)
				{
					int val = GetData2(C, k, j, dims[0], dims[1], (int)0);
					if (val > 0 && count[val] < minsize)
					{
						SetData4(F, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0);
						nremoved++;
					}
				}
			}
		}
	}
	printf("removeIsolated: %d voxels removed.\n", nremoved);
}

void
ExtactSurface(vector<unsigned char>& T, 
			  const vector<int>& C, 
			  const vector<unsigned char>& L, 
			  const int* dims)
{
	for (int m = 0; m<dims[3]; ++m)
	{
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				for (int k = 0; k < dims[0]; ++k)
				{
					if (onSurface4(L, k, j, i, m, dims))
					{
						int cval = GetData4(C, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (int)0);
						if (cval)
						{
							SetData4(T, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1);
						}
					}
				}
			}
		}
	}
}

vector<CParticle4D>
collectSphericalFront(vector<unsigned char>& F,
					  const vector<unsigned char>& S, //foreground
					  const vector<float>& D,
					  const vector<CParticle4D>& vp, 
					  const vector<float>& v,
					  const int* dims)
{
	vector<unsigned char> L(S.size(), 0);

	for (int m = 0; m<dims[3]; ++m)
	{
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				for (int k = 0; k < dims[0]; ++k)
				{
					for (int n = 0; n < vp.size(); ++n)
					{
						int x = vp[n].m_X;
						int y = vp[n].m_Y;
						int z = vp[n].m_Z;
						int t = vp[n].m_T;
						float dval = GetData4(D, x, y, z, t, dims[0], dims[1], dims[2], dims[3], (float)0);
						float dx = (k - x)*v[0];
						float dy = (j - y)*v[1];
						float dz = (i - z)*v[2];
						float dt = (m - t)*v[3];
						if (dx*dx + dy*dy + dz*dz + dt*dt < dval*dval)
						{
							SetData4(L, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1);
						}
					}
				}
			}
		}
	}
	vector<CParticle4D> vfront;
	for (int m = 0; m<dims[3]; ++m)
	{
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				for (int k = 0; k < dims[0]; ++k)
				{
					if (GetData4(L, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0))
					{
						int sval = GetData4(S, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0);
						if (sval)
						{
							SetData4(F, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1);
							if (onSurface4(L, k, j, i, m, dims))
							{
								vfront.push_back(CParticle4D(k, j, i, m));
							}
						}
					}
				}
			}
		}
	}
	return vfront;
}

int
FreeGrowthHere(vector<int>& B,
		   const vector<float>& D,
		   const vector<CParticle4D>& vSeeds,
		   bool bStrict,
		   const int* dims)
{
	set<CParticle4D> vp;
	for (int i = 0; i < vSeeds.size(); ++i)
	{
		vp.insert(vSeeds[i]);
	}
	int cnt = 1;
	while(!vp.empty())
	{
		set<CParticle4D> vp2;
		for (set<CParticle4D>::iterator it = vp.begin(); it != vp.end(); ++it)
		{
			SetVoxel(B, *it, (int)cnt, dims);
		}
		for (set<CParticle4D>::iterator it = vp.begin(); it != vp.end(); ++it)
		{
			CParticle4D p4 = *it;
			int x1 = p4.m_X;
			int y1 = p4.m_Y;
			int z1 = p4.m_Z;
			int t1 = p4.m_T;
			float dval1 = GetData4(D, x1, y1, z1, t1, dims[0], dims[1], dims[2], dims[3], (float)0);
			for(int j=0; j<NumNeighborsHere; ++j)
			{
				int x2 = x1 + XOffsetHere[j];
				int y2 = y1 + YOffsetHere[j];
				int z2 = z1 + ZOffsetHere[j];
				int t2 = t1 + TOffsetHere[j];
				if (!GetData4(B, x2, y2, z2, t2, dims[0], dims[1], dims[2], dims[3], (int)1))
				{
					float dval2 = GetData4(D, x2, y2, z2, t2, dims[0], dims[1], dims[2], dims[3], (float)0);
					if(dval2>0)
					{
						if((bStrict && dval1 > dval2) || (!bStrict && dval1 >= dval2))
						{
							CParticle4D pn(x2, y2, z2, t2);
							vp2.insert(pn);
						}
					}
				}
			}
		}
		vp = vp2;
		cnt++;
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
	for (int m = 0; m < dims[3]; ++m)
	{
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				for (int k = 0; k < dims[0]; ++k)
				{
					unsigned char tval = GetData4(T, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0);
					if (tval)
					{
						int sval = GetData4(S, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (int)0);
						/*if (sval == 1)
						{
							//this is within the initial sphere.  Accept it.
							SetData4(U, k, j, i, m, dims[0], dims[1], dims[2], dims[3], tval);
							continue;
						}*/
						float dval = GetData4(D, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (float)0);
						float maxrate = 0;
						vector<CParticle4D> vsource;
						for (int n = 0; n<NumNeighborsHere; ++n)
						{
							int x2 = k + XOffsetHere[n];
							int y2 = j + YOffsetHere[n];
							int z2 = i + ZOffsetHere[n];
							int t2 = m + TOffsetHere[n];
							float dval2 = GetData4(D, x2, y2, z2, t2, dims[0], dims[1], dims[2], dims[3], (float)0);
							float len = sqrt((float)XOffsetHere[m] * XOffsetHere[m] + YOffsetHere[m] * YOffsetHere[m] + ZOffsetHere[m] * ZOffsetHere[m] + TOffsetHere[m] * TOffsetHere[m]);
							float rate = (dval2 - dval) / len;
							if (rate > maxrate)
							{
								maxrate = rate;
								vsource.clear();
								vsource.push_back(CParticle4D(x2, y2, z2, t2));
							}
							else if (rate == maxrate)
							{
								vsource.push_back(CParticle4D(x2, y2, z2, t2));
							}
						}
						bool bSuccess = false;
						//if (maxrate > 0.75)
						{
							for (int n = 0; n < vsource.size(); ++n)
							{
								int sval2 = GetVoxel(S, vsource[n], (int)0, dims);
								if (sval2 == sval - 1)
								{
									bSuccess = true;
									break;
								}
							}
						}
						if (bSuccess)
						{
							SetData4(U, k, j, i, m, dims[0], dims[1], dims[2], dims[3], tval);
						}
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
	vector<int> nbh = MakeEightNeighborhood(4, dims); //MakeEightNeighborhood(4, dims);
	int nc = ConnectedComponentAnalysisBigger(C, T, nbh, (unsigned char)0, 4, dims); 
	printf("There are %d clusters in T.\n", nc);
	int minLabel = std::numeric_limits<int>::max();
	for (int i = 0; i < S.size(); ++i)
	{
		if (S[i] && U[i])
		{
			minLabel = Min(minLabel, S[i]);
		}
	}
	printf("minlabel = %d.\n", minLabel);
	for(int n=1; n<=nc; ++n)
	{
		vector<unsigned char> L(U.size(), 0);
		for(int i=0; i<L.size(); ++i)
		{
			if(C[i]==n && U[i]) ///Including U[i] correct?
			{
				L[i]=1;
			}
		}
		vector<int> D(S.size(), 0);
		int nc2 = ConnectedComponentAnalysisBigger(D, L, nbh, (unsigned char)0, 4, dims); 
		printf("minlabel = %d.  There are %d clusters in Component %d.\n", minLabel, nc2, n);

		vector<bool> bFlag(nc2 + 1, false);
		vector<int> counts(nc2 + 1, 0);
		for (int i = 0; i<D.size(); ++i)
		{
			if(D[i]>0)
			{
				if(S[i] == minLabel)
				{
					bFlag[D[i]]=true;
				}
			}
			counts[D[i]]++;
		}
		for(int i=1; i<bFlag.size(); ++i)
		{
			/*if(bFlag[i])
			{
				printf("Component %d is accepted.\n", i);
			}*/
			//printf("%d %s %d\n", i, bFlag[i] ? "True" : "False", counts[i]);
		}
		//printf("\n");
		for(int i=0; i<D.size(); ++i)
		{
			if(D[i]>0 && bFlag[D[i]])
			{
				V[i] = 1;
			}
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 3 || nlhs < 0)
	{
		//mexErrMsgTxt("Usage: [C [S D]] = TraceTest(S, L, seed)");
		mexErrMsgTxt("Usage: [C] = TraceTest(L, v)");
		return;
	}

	initializeOffsetArrays(); //need to initialize those...

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

	vector<CParticle4D> vpclick; //(1, CParticle(vseed[0], vseed[1], vseed[2]));
	for(int i=0; i<vseed.size(); i+=4)
	{
		printf("straightGrowth: Seed = (%d, %d, %d, %d)\n", vseed[i], vseed[i+1], vseed[i+2], vseed[i+3]);
		vpclick.push_back(CParticle4D(vseed[i], vseed[i + 1], vseed[i + 2], vseed[i + 3]));
	}

	vector<float> D(nvoxels,0);
	vector<unsigned char> iL(nvoxels,0);
	for(int i=0; i<nvoxels; ++i)
	{
		iL[i] = L[i] ? 0: 1;
	}
	DistanceTransformEuclidF(D, iL, v, ndimL, dimsL);

	vector<unsigned char> F(nvoxels, 0);
	vector<CParticle4D> vfront = collectSphericalFront(F, L, D, vpclick, v, dimsL);

	printf("There are %d front voxels.\n", vfront.size());
	vector<int> S(nvoxels, 0);
	for(int i=0; i<nvoxels; ++i)
	{
		if(F[i]) S[i] = 1;
	}
	FreeGrowthHere(S, D, vfront, true, dimsL);		

	vector<unsigned char> T(nvoxels, 0);
	ExtactSurface(T, S, L, dimsL);

	vector<unsigned char> T2(T.size(), 0);
	CheckTheGrowth(T2, T, S, D, dimsL);

	removeIsolated(T2, 3, dimsL);

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

