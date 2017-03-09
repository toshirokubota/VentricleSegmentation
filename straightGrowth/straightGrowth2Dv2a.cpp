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
#include <szConvolutionTemplate.h>
//#include <szFigureGroundSeparation.h>
//#include <szExtractForeground.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>
//#include <szConvexHull3D.h>
//#include <szCutOutNodule.h>
#include <szParticle4D.h>
#include <szDistanceTransform.h>
#include <szConnectedComponent.h>
#include <szLocalExtrema.h>
//#include <szTrace.h>
//#include <szRetrace.h>
#include <szDistanceTransformNonIsotropic.h>
//#define FOUR_NEIGHBORHOOD
#include <Graph.h>
#include <GraphFactory.h>
#include <Kruskal.h>

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
const int NumNeighborsHere = 8; // 3 * NumNeighbors + 2;
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
	for (int j = -1; j <= 1; ++j)
	{
		for (int k = -1; k <= 1; ++k)
		{
			if (j == 0 && k == 0)
			{
				continue;
			}
			else
			{
				XOffsetArray[idx] = k;
				YOffsetArray[idx] = j;
				OffsetLength[idx] = sqrt(k*k + j*j);
				idx++;
			}
		}
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
	return SetData2(A, p.m_X, p.m_Y, dims[0], dims[1], value);
}

template<class T>
T
GetVoxel(const vector<T>& A,
const CParticle4D& p,
T defaultValue,
const int* dims)
{
	return GetData2(A, p.m_X, p.m_Y, dims[0], dims[1], defaultValue);
}

/*
onSurface3 with a general neighborhood
*/
template<class T>
bool
onSurface2(const vector<T>& A, int x, int y, const int* dims) 
{
	if (GetData2(A, x, y, dims[0], dims[1], (T)0))
	{
		for(int m=0; m<NumNeighborsHere; ++m)
		{
			int x2 = x + XOffsetHere[m];
			int y2 = y + YOffsetHere[m];
			if (GetData2(A, x2, y2, dims[0], dims[1], (T)1) == 0)
			{
				return true;
			}
		}
	}
	return false;
}

vector<CParticle4D>
ExtactSurface(const vector<unsigned char>& L, 
			  const int* dims)
{
	vector<CParticle4D> front;
	for (int j = 0; j < dims[1]; ++j)
	{
		for (int k = 0; k < dims[0]; ++k)
		{
			if (onSurface2(L, k, j, dims))
			{
				front.push_back(CParticle4D(k, j));
			}
		}
	}
	return front;
}

float 
Length2(float x, float y)
{
	return sqrt(x*x + y*y);
}

float*
Gradient2(const vector<float>& D,
const CParticle4D& p,
const int* dims)
{
	int x = p.m_X;
	int y = p.m_Y;
	float dval1 = GetData2(D, x, y, dims[0], dims[1], (float)0);
	float dvalW = GetData2(D, x - 1, y, dims[0], dims[1], (float)dval1);
	float dvalE = GetData2(D, x + 1, y, dims[0], dims[1], (float)dval1);
	float dvalN = GetData2(D, x, y - 1, dims[0], dims[1], (float)dval1);
	float dvalS = GetData2(D, x, y + 1, dims[0], dims[1], (float)dval1);
	float* g = new float[2];
	g[0] = dvalE - dvalW;
	g[1] = dvalS - dvalN;
	float len = Length2(g[0], g[1]);
	g[0] /= len; g[1] /= len;
	return g;
}

map<CParticle4D,CParticle4D>
SteepestGrowthHere(vector<CParticle4D>& front,
		   const vector<float>& D,
		   const int* dims)
{
	int label = 1;
	map<CParticle4D, CParticle4D> coremap;
	map<CParticle4D, bool> bmap;
	map<CParticle4D, float*> grad;
	for (vector<CParticle4D>::iterator it = front.begin(); it != front.end(); ++it)
	{
		coremap[*it] = *it;
		bmap[*it] = false;
		grad[*it] = Gradient2(D, *it, dims);
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
			float dval1 = GetVoxel(D, p4, 0.0f, dims);
			float* g = grad[p0];
			float maxd = 0; // g[0] * (x1 - p0.m_X) + g[1] * (y1 - p0.m_Y);
			CParticle4D cand = p4;
			bool bfound = false;
			if (p0.m_X + 1 == 55 && p0.m_Y + 1 == 37)
			{
				dval1 += 0;
			}
			for (int j = 0; j<NumNeighborsHere; ++j)
			{
				int x2 = x1 + XOffsetHere[j];
				int y2 = y1 + YOffsetHere[j];
				float d = (g[0] * (x2 - p0.m_X) + g[1] * (y2 - p0.m_Y)) / OffsetLength[j];
				if (d > maxd)
				{
					maxd = d;
					cand = CParticle4D(x2, y2);
				}
			}
			float dval2 = GetVoxel(D, cand, 0.0f, dims);
			if (dval2 > dval1)
			{
				coremap[p0] = cand;
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
	for (vector<CParticle4D>::iterator it = front.begin(); it != front.end(); ++it)
	{
		delete grad[*it];
	}
	return coremap;
}

bool
isNeighbor(const CParticle4D& a, const CParticle4D& b, int space=1)
{
	/*if (a.m_X == b.m_X && a.m_Y == b.m_Y) return false;
	else if (Abs(a.m_X - b.m_X) <= 1 && a.m_Y == b.m_Y) return true;
	else if (Abs(a.m_Y - b.m_Y) <= 1 && a.m_X == b.m_X) return true;
	else return false;*/
	if (a.m_X == b.m_X && a.m_Y == b.m_Y) return false;
	else
	{
		for (int y = a.m_Y - space; y <= a.m_Y + space; ++y)
		{
			for (int x = a.m_X - space; x <= a.m_X + space; ++x)
			{
				if (y == b.m_Y && x == b.m_X)
				{
					return true;
				}

			}
		}
		return false;
	}
}

vector<Vertex<CParticle4D>*>
makeGraph(set<CParticle4D>& core,
map<CParticle4D, CParticle4D>& pmap)
{
	GraphFactory<CParticle4D>& factory = GraphFactory<CParticle4D>::GetInstance();
	vector<Vertex<CParticle4D>*> vertices;

	map<CParticle4D, set<CParticle4D>> revmap;
	for (set<CParticle4D>::iterator it = core.begin(); it != core.end(); ++it) {
		revmap[*it] = set<CParticle4D>();
		vertices.push_back(factory.makeVertex(*it));
	}
	for (map<CParticle4D, CParticle4D>::iterator it = pmap.begin(); it != pmap.end(); ++it) {
		revmap[it->second].insert(it->first);
	}
	for (int i = 0; i < vertices.size(); ++i)
	{
		Vertex<CParticle4D>* u = vertices[i];
		if (u->key.m_X == 60 && u->key.m_Y == 46)
		{
			i += 0;
		}
		for (int j = i + 1; j < vertices.size(); ++j)
		{
			Vertex<CParticle4D>* v = vertices[j];
			bool bfound = false;
			for (set<CParticle4D>::iterator it = revmap[u->key].begin(); it != revmap[u->key].end() && bfound==false; ++it)
			{
				for (set<CParticle4D>::iterator it2 = revmap[v->key].begin(); it2 != revmap[v->key].end() && bfound == false; ++it2)
				{
					if (isNeighbor(*it, *it2))
					{
						bfound = true;
					}
				}
			}
			if (bfound)
			{
				float w = Length2(u->key.m_X - v->key.m_X, u->key.m_Y - v->key.m_Y);
				Edge<CParticle4D>* e1 = factory.makeEdge(u, v, w);
				Edge<CParticle4D>* e2 = factory.makeEdge(v, u, w);
				u->Add(e1);
				v->Add(e2);
			}
		}
	}

	return vertices;
}

vector<int>
clusterVertices(vector<Vertex<CParticle4D>*>& vertices, float thres)
{
	vector<Node<int>*> nodes;
	map<Vertex<CParticle4D>*,int> imap;
	for (int i = 0; i < vertices.size(); ++i)
	{
		nodes.push_back(makeset(i));
		imap[vertices[i]] = i;
	}
	for (int i = 0; i < vertices.size(); ++i)
	{
		for (int j = 0; j < vertices[i]->aList.size(); ++j)
		{
			Edge<CParticle4D>* ed = vertices[i]->aList[j];
			if (ed->w < thres)
			{
				int k = imap[ed->v];
				merge(nodes[i], nodes[k]);
			}
		}
	}
	vector<Node<int>*> reps = clusters(nodes);
	printf("clusterCore: found %d clusters.\n", reps.size());
	vector<int> label(nodes.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = nodes[i]->key;
		int m = distance(reps.begin(), find(reps.begin(), reps.end(), findset(nodes[i])));
		label[k] = m + 1; //make it positive
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return label;
}

map<CParticle4D,int>
partitionSurface(vector<vector<CParticle4D>>& cores,
map<CParticle4D, CParticle4D>& coremap)
{
	map<CParticle4D, int> lmap;
	for (int i = 0; i < cores.size(); ++i) {
		for (int j = 0; j < cores[i].size(); ++j)
		{
			lmap[cores[i][j]] = i;
		}
	}
	map<CParticle4D, int> label;
	for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it)
	{
		label[it->first] = lmap[it->second];
	}
	return label;
}

void
partitionRegion(vector<int>& S, //labeled image
				const vector<unsigned char>& L, //foreground
				map<CParticle4D, int>& partition, 
				int minsize,
				int ndim, const int* dims)
{
	int numclusters = 0;
	for (map<CParticle4D, int>::iterator it = partition.begin(); it != partition.end(); ++it)
	{
		numclusters = Max(numclusters, it->second + 1);
	}
	vector<vector<CParticle4D>> Q(numclusters);
	for (map<CParticle4D, int>::iterator it = partition.begin(); it != partition.end(); ++it)
	{
		Q[it->second].push_back(it->first);
	}
	for (int i = 0; i < Q.size(); ++i)
	{
		if (Q[i].size() < minsize)
		{
			Q[i].clear();
		}
	}
	while (true)
	{
		bool bcontinue = false;
		map<CParticle4D, int> marks;
		for (int i = 0; i < Q.size(); ++i)
		{
			for (int j = 0; j < Q[i].size(); ++j)
			{
				CParticle4D p = Q[i][j];
				for (int n = 0; n < NumNeighborsHere; ++n)
				{
					CParticle4D q(p.m_X + XOffsetHere[n], p.m_Y + YOffsetHere[n]);
					if (GetVoxel(L, q, (unsigned char)0, dims))
					{
						if (GetVoxel(S, q, 0, dims) == 0)
						{
							if (marks.find(q) == marks.end())
							{
								marks[q] = i;
							}
							else
							{
								if (marks[q] != i)
								{
									marks[q] = -1; //conflict
								}
							}
							bcontinue = true;
						}
					}
				}
			}
		}
		vector<vector<CParticle4D>> Q2(numclusters);
		for (map<CParticle4D, int>::iterator it = marks.begin(); it != marks.end(); ++it)
		{
			int m = it->second;
			SetVoxel(S, it->first, m, dims);
			if (m >= 0)
			{
				Q2[m].push_back(it->first);
			}
		}
		Q = Q2;

		if (bcontinue == false) break;
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1|| nlhs < 0)
	{
		//mexErrMsgTxt("Usage: [C [S D]] = TraceTest(S, L, seed)");
		mexErrMsgTxt("Usage: [C] = straightGrowth(L)");
		return;
	}

	initializeOffsetArrays(); //need to initialize those...

	//load figure-ground segmentation
	int ndimL;
	const int* dimsL;
	mxClassID classL;
	vector<unsigned char> L;
	LoadData(L,prhs[0],classL,ndimL,&dimsL);

	float thres = 2.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[1], classMode);
	}
	int minsize = 0.0f;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minsize, prhs[2], classMode);
	}
	printf("threshold = %f, minsize = %d\n", thres, minsize);

	int nvoxels = numberOfElements(ndimL,dimsL);

	vector<CParticle4D> vfront = ExtactSurface(L, dimsL);
	printf("There are %d boundary voxels.\n", vfront.size());

	GraphFactory<CParticle4D>& factory = GraphFactory<CParticle4D>::GetInstance();
	vector<float> D0(nvoxels, 0);
	vector<unsigned char> iL(nvoxels,0);
	for(int i=0; i<nvoxels; ++i)
	{
		iL[i] = L[i] ? 0: 1;
	}
	vector<float> v(2, 1.0f);
	DistanceTransformEuclidF(D0, iL, v, ndimL, dimsL);
	vector<float> filter = makeGaussianFilter(1.0f, 3.0f, 9);
	vector<float> D(nvoxels, 0);
	separableConvolution(D, D0, filter, CBE_ZERO, ndimL, dimsL);

	map<CParticle4D,CParticle4D> coremap = SteepestGrowthHere(vfront, D, dimsL);

	set<CParticle4D> coreset;
	for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it)
	{
		coreset.insert(it->second);
	}
	
	vector<Vertex<CParticle4D>*> vertices = makeGraph(coreset, coremap);
	vector<int> label = clusterVertices(vertices, thres);

	int numclusters = 0;
	for (int i = 0; i < label.size(); ++i)
	{
		numclusters = Max(numclusters, label[i] + 1);
	}
	vector<vector<CParticle4D>> cores(numclusters);
	for (int i = 0; i < vertices.size(); ++i)
	{
		cores[label[i]].push_back(vertices[i]->key);
	}

	map<CParticle4D, int> partition = partitionSurface(cores, coremap);

	vector<int> S(nvoxels, 0);
	partitionRegion(S, L, partition, minsize, ndimL, dimsL);

	if(nlhs >= 1)
	{
		plhs[0] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 2)
	{
		int dims[] = { coremap.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it, i++)
		{
			SetData2(F, i, 0, dims[0], dims[1], it->first.m_X);
			SetData2(F, i, 1, dims[0], dims[1], it->first.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], it->second.m_X);
			SetData2(F, i, 3, dims[0], dims[1], it->second.m_Y);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		int dims[] = { factory.edges.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], (float)factory.edges[i]->u->key.m_X);
			SetData2(F, i, 1, dims[0], dims[1], (float)factory.edges[i]->u->key.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], (float)factory.edges[i]->v->key.m_X);
			SetData2(F, i, 3, dims[0], dims[1], (float)factory.edges[i]->v->key.m_Y);
			SetData2(F, i, 4, dims[0], dims[1], (float)factory.edges[i]->w);
		}
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(D, mxSINGLE_CLASS, ndimL, dimsL);
	}
	factory.Clean();
	mexUnlock();
}

