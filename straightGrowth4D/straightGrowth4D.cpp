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
#include <DisjointSet.h>
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
float OffsetLength[NumNeighborsHere];

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
						OffsetLength[idx] = sqrt(k*k + j*j + i*i + m*m);
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

set<CParticle4D>
FreeGrowthHere(vector<int>& B,
const vector<float>& D,
const vector<CParticle4D>& vSeeds,
bool bStrict,
const int* dims)
{
	set<CParticle4D> front;
	set<CParticle4D> vp;
	for (int i = 0; i < vSeeds.size(); ++i)
	{
		vp.insert(vSeeds[i]);
	}
	int cnt = 1;
	while (!vp.empty())
	{
		set<CParticle4D> vp2;
		for (set<CParticle4D>::iterator it = vp.begin(); it != vp.end(); ++it)
		{
			SetVoxel(B, *it, (int)cnt, dims);
			if (Abs(GetVoxel(D, *it, 0.0f, dims) - 1.0f) < 0.01)
			{
				front.insert(*it);
			}
		}
		for (set<CParticle4D>::iterator it = vp.begin(); it != vp.end(); ++it)
		{
			CParticle4D p4 = *it;
			int x1 = p4.m_X;
			int y1 = p4.m_Y;
			int z1 = p4.m_Z;
			int t1 = p4.m_T;
			float dval1 = GetData4(D, x1, y1, z1, t1, dims[0], dims[1], dims[2], dims[3], (float)0);
			for (int j = 0; j<NumNeighborsHere; ++j)
			{
				int x2 = x1 + XOffsetHere[j];
				int y2 = y1 + YOffsetHere[j];
				int z2 = z1 + ZOffsetHere[j];
				int t2 = t1 + TOffsetHere[j];
				if (!GetData4(B, x2, y2, z2, t2, dims[0], dims[1], dims[2], dims[3], (int)1))
				{
					float dval2 = GetData4(D, x2, y2, z2, t2, dims[0], dims[1], dims[2], dims[3], (float)0);
					if (dval2>0)
					{
						if ((bStrict && dval1 > dval2) || (!bStrict && dval1 >= dval2))
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
	printf("FreeGrowth: There are %d front voxels.\n", front.size());
	return front;
}

float 
Length4(float x, float y, float z, float t)
{
	return sqrt(x*x + y*y + z*z + t*t);
}

float*
Gradient4(const vector<float>& D,
const CParticle4D& p,
const int* dims)
{
	int x = p.m_X;
	int y = p.m_Y;
	int z = p.m_Z;
	int t = p.m_T;
	float dval1 = GetData4(D, x, y, z, t, dims[0], dims[1], dims[2], dims[3], (float)0);
	float dvalW = GetData4(D, x - 1, y, z, t, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float dvalE = GetData4(D, x + 1, y, z, t, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float dvalN = GetData4(D, x, y - 1, z, t, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float dvalS = GetData4(D, x, y + 1, z, t, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float dvalB = GetData4(D, x, y, z - 1, t, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float dvalF = GetData4(D, x, y, z + 1, t, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float dvalPrev = GetData4(D, x, y, z, t - 1, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float dvalNext = GetData4(D, x, y, z, t + 1, dims[0], dims[1], dims[2], dims[3], (float)dval1);
	float* g = new float[4];
	g[0] = dvalE - dvalW;
	g[1] = dvalS - dvalN;
	g[2] = dvalF - dvalB;
	g[3] = dvalNext - dvalPrev;
	float len = Length4(g[0], g[1], g[2], g[3]);
	g[0] /= len; g[1] /= len; g[2] /= len; g[3] /= len;
	return g;
}

float projectedDistanceValue(const vector<float>& D, CParticle4D& p0, CParticle4D& p, float* grd, const int* dims)
{
	return 0.0f;
}

map<CParticle4D,CParticle4D>
SteepestGrowthHere(vector<int>& B,
		   const vector<float>& D,
		   const set<CParticle4D>& front,
		   const int* dims)
{
	int label = 1;
	map<CParticle4D, CParticle4D> coremap;
	map<CParticle4D, bool> bmap;
	map<CParticle4D, float*> grad;
	for (set<CParticle4D>::iterator it = front.begin(); it != front.end(); ++it)
	{
		coremap[*it] = *it;
		bmap[*it] = false;
		SetVoxel(B, *it, label, dims);
		grad[*it] = Gradient4(D, *it, dims);
	}
	int total = 0;
	bool bchanged = true;
	while (bchanged)
	{
		bchanged = false;
		for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it)
		{
			CParticle4D p0 = it->first;
			if (bmap[p0] == true) continue;
			float dval0 = GetVoxel(D, p0, 0.0f, dims);

			CParticle4D p4 = it->second;
			int x1 = p4.m_X;
			int y1 = p4.m_Y;
			int z1 = p4.m_Z;
			int t1 = p4.m_T;
			float dval1 = GetVoxel(D, p4, 0.0f, dims);
			float* g = grad[p0];
			float maxd = g[0] * (x1 - p0.m_X) + g[1] * (y1 - p0.m_Y) + g[2] * (z1 - p0.m_Z) + g[3] * (t1 - p0.m_T);
			CParticle4D cand = p4;
			bool bfound = false;
			for (int j = 0; j<NumNeighborsHere; ++j)
			{
				int x2 = x1 + XOffsetHere[j];
				int y2 = y1 + YOffsetHere[j];
				int z2 = z1 + ZOffsetHere[j];
				int t2 = t1 + TOffsetHere[j];
				float d = g[0] * (x2 - p0.m_X) + g[1] * (y2 - p0.m_Y) + g[2] * (z2 - p0.m_Z) + g[3] * (t2 - p0.m_T);
				if (d > maxd)
				{
					maxd = d;
					cand = CParticle4D(x2, y2, z2, t2);
				}
			}
			float dval2 = GetVoxel(D, cand, 0.0f, dims);
			if (dval2 > dval1)
			{
				coremap[p0] = cand;
				SetVoxel(B, cand, label, dims);
				bchanged = true;
			}
			else
			{
				bmap[p0] = true;
			}
		}
		label++;
	}
	printf("Done tracing %d voxels.\n", coremap.size());
	for (set<CParticle4D>::iterator it = front.begin(); it != front.end(); ++it)
	{
		delete grad[*it];
	}
	return coremap;
}

map<CParticle4D,unsigned int>
clusterCores(set<CParticle4D>& core,
const vector<float>& D,
float proportion,
const int* dims)
{
	vector<CParticle4D> vcore;
	vector<Node<unsigned int>*> nodes;
	for (set<CParticle4D>::iterator it = core.begin(); it != core.end(); ++it)
	{
		vcore.push_back(*it);
		nodes.push_back(makeset(vcore.size() - 1));
	}
	for (int i = 0; i < vcore.size(); ++i)
	{
		CParticle4D p1 = vcore[i];
		float dval1 = GetVoxel(D, p1, 0.0f, dims);
		for (int j = i + 1; j < vcore.size(); ++j)
		{
			CParticle4D p2 = vcore[j];
			float dval2 = GetVoxel(D, p2, 0.0f, dims);
			float sep = sqrt((p1.m_X - p2.m_X)*(p1.m_X - p2.m_X) + (p1.m_Y - p2.m_Y) * (p1.m_Y - p2.m_Y) + (p1.m_Z - p2.m_Z)*(p1.m_Z - p2.m_Z) + (p1.m_T - p2.m_T)*(p1.m_T - p2.m_T));
			bool bmerge = false;
			if (sep < proportion * sqrt(dval1 * dval2))
			{
				merge(nodes[i], nodes[j]);
				bmerge = true;
			}
			if (p2.m_X + 1 == 58 && p2.m_Y + 1 == 50 && p2.m_Z + 1 == 52 && p2.m_T + 1 == 17)
			{
				if (p1.m_X + 1 == 59 && p1.m_Y + 1 == 49 && p1.m_Z + 1 == 52 && p1.m_T + 1 == 17)
				{
					printf("dval1=%f, dval2=%f, sep=%f, merge=%d\n", dval1, dval2, sep, bmerge);
				}
			}
		}
	}
	vector<Node<unsigned int>*> reps = clusters(nodes);
	printf("clusterCore: found %d clusters.\n", reps.size());
	map<CParticle4D, unsigned int> pmap;
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = nodes[i]->key;
		int m = distance(reps.begin(), find(reps.begin(), reps.end(), findset(nodes[i])));
		pmap[vcore[k]] = m;
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return pmap;
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

	float proportion = 0.5;
	if (nrhs >= 4)
	{
		mxClassID classMode;
		ReadScalar(proportion, prhs[3], classMode);
	}
	printf("proportion = %f\n", proportion);

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

	printf("There are %d spherical front voxels.\n", vfront.size());
	vector<int> S(nvoxels, 0);
	for(int i=0; i<nvoxels; ++i)
	{
		if(F[i]) S[i] = 1;
	}
	set<CParticle4D> surf = FreeGrowthHere(S, D, vfront, true, dimsL);		
	vector<int> S2(nvoxels, 0);
	map<CParticle4D,CParticle4D> coremap = SteepestGrowthHere(S2, D, surf, dimsL);

	set<CParticle4D> coreset;
	for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it)
	{
		coreset.insert(it->second);
	}
	map<CParticle4D, unsigned int> labels = clusterCores(coreset, D, proportion, dimsL);

	//select cores that are inside the initial spherical area
	int numclusters = 0;
	for (map<CParticle4D, unsigned int>::iterator it = labels.begin(); it != labels.end(); ++it)
	{
		numclusters = Max(it->second + 1, numclusters);
	}
	vector<bool> bSelection(numclusters, false);
	for (set<CParticle4D>::iterator it = coreset.begin(); it != coreset.end(); ++it)
	{
		if (GetVoxel(F, *it, (unsigned char)0, dimsL))
		{
			int m = labels[*it];
			bSelection[m] = true;
		}
	}

	vector<unsigned char> T(nvoxels, 0);
	vector<unsigned char> J(nvoxels, 0);
	for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it)
	{
		int m = labels[it->second];
		SetVoxel(J, it->first, (unsigned char)1, dimsL);
		if (bSelection[m])
		{
			SetVoxel(T, it->first, (unsigned char)1, dimsL);
			SetVoxel(J, it->second, (unsigned char)3, dimsL);
		}
		else
		{
			SetVoxel(J, it->second, (unsigned char)2, dimsL);
		}
	}

	if(nlhs >= 1)
	{
		plhs[0] = StoreData(S,mxINT32_CLASS,ndimL,dimsL);
	}
	if(nlhs >= 2)
	{
		plhs[1] = StoreData(T,mxUINT8_CLASS,ndimL,dimsL);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(S2, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(D, mxSINGLE_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 5)
	{
		plhs[4] = StoreData(J, mxUINT8_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 6)
	{
		int dims[] = { coremap.size(), 8 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0; 
		for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it, i++)
		{
			SetData2(F, i, 0, dims[0], dims[1], it->first.m_X);
			SetData2(F, i, 1, dims[0], dims[1], it->first.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], it->first.m_Z);
			SetData2(F, i, 3, dims[0], dims[1], it->first.m_T);
			SetData2(F, i, 4, dims[0], dims[1], it->second.m_X);
			SetData2(F, i, 5, dims[0], dims[1], it->second.m_Y);
			SetData2(F, i, 6, dims[0], dims[1], it->second.m_Z);
			SetData2(F, i, 7, dims[0], dims[1], it->second.m_T);
		}
		plhs[5] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
}

