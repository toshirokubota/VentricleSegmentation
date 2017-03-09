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

vector<float>
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
	vector<float> g(2);
	g[0] = dvalE - dvalW;
	g[1] = dvalS - dvalN;
	float len = Length2(g[0], g[1]);
	g[0] /= len; g[1] /= len;
	return g;
}

/*
Trace a path along the gradient (in D) at the p, until it can no longer ascend (in D) 
or one of target pixels is reached.
*/
CParticle4D
steepestAscent(CParticle4D p, set<CParticle4D>& target, const vector<float>& D, const int* dims)
{
	vector<float> g = Gradient2(D, p, dims);
	CParticle4D p0 = p;
	while (true)
	{
		CParticle4D cand;
		float dval1 = GetVoxel(D, p, 0.0f, dims);
		float maxd = 0.0f;
		for (int j = 0; j<NumNeighborsHere; ++j)
		{
			int x2 = p.m_X + XOffsetHere[j];
			int y2 = p.m_Y + YOffsetHere[j];
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
			p = cand;
			if (target.find(p) != target.end())
			{
				break; //target is reached.
			}
		}
		else
		{
			break;
		}
	}
	return p;
}

/*
From each surface pixel (FRONT), find core pixels by going up straight in the distance map (D).
*/
map<CParticle4D, CParticle4D>
SteepestGrowthHere(vector<CParticle4D>& front,
		   const vector<float>& D,
		   const int* dims)
{
	map<CParticle4D, CParticle4D> coremap;
	set<CParticle4D> target; //dummy - empty set.
	for (vector<CParticle4D>::iterator it = front.begin(); it != front.end(); ++it)
	{
		coremap[*it] = steepestAscent(*it, target, D, dims);
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


/*
Connect cores that are within THRES distance away.
*/
vector<Vertex<CParticle4D>*>
makeGraph(map<CParticle4D,CParticle4D>& coremap, float thres)
{
	GraphFactory<CParticle4D>& factory = GraphFactory<CParticle4D>::GetInstance();
	set<CParticle4D> core;
	for (map<CParticle4D, CParticle4D>::iterator it = coremap.begin(); it != coremap.end(); ++it)
	{
		core.insert(it->second);
	}
	vector<Vertex<CParticle4D>*> vertices;

	for (set<CParticle4D>::iterator it = core.begin(); it != core.end(); ++it) {
		vertices.push_back(factory.makeVertex(*it));
	}
	for (int i = 0; i < vertices.size(); ++i)
	{
		Vertex<CParticle4D>* u = vertices[i];
		for (int j = i + 1; j < vertices.size(); ++j)
		{
			Vertex<CParticle4D>* v = vertices[j];
			float w = Length2(u->key.m_X - v->key.m_X, u->key.m_Y - v->key.m_Y);
			if (w < thres)
			{
				Edge<CParticle4D>* e1 = factory.makeEdge(u, v, w);
				Edge<CParticle4D>* e2 = factory.makeEdge(v, u, w);
				u->Add(e1);
				v->Add(e2);
			}
		}
	}

	return vertices;
}

/*
Cluster vertices by a simple connected component.
*/
map<CParticle4D,int>
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
			int k = imap[ed->v];
			merge(nodes[i], nodes[k]);
		}
	}
	vector<Node<int>*> reps = clusters(nodes);
	printf("clusterCore: found %d clusters.\n", reps.size());
	map<CParticle4D,int> lmap;
	for (int i = 0; i < nodes.size(); ++i)
	{
		int m = distance(reps.begin(), find(reps.begin(), reps.end(), findset(nodes[i])));
		lmap[vertices[i]->key] = m + 1; //make it positive.
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return lmap;
}

map<CParticle4D,int>
partitionSurface(vector<CParticle4D>& surface,
map<CParticle4D,int>& label,
vector<float>& D, const int* dims)
{
	set<CParticle4D> coreset;
	for (map<CParticle4D, int>::iterator it = label.begin(); it != label.end(); ++it)
	{
		coreset.insert(it->first);
	}
	map<CParticle4D, int> lmap;
	for (int i = 0; i < surface.size(); ++i)
	{
		CParticle4D p = surface[i];
		CParticle4D core = steepestAscent(p, coreset, D, dims);
		if (label.find(core) != label.end())
		{
			lmap[p] = label[core];
		}
		else
		{
			lmap[p] = -1; //this should not happen...
		}
	}
	return lmap;
}

void
partitionRegion(vector<int>& S, //labeled image
				const vector<unsigned char>& L, //foreground
				const vector<float>& D,
				map<CParticle4D, int>& cores, 
				int minsize,
				int ndim, const int* dims)
{
	int numclusters = 0;
	for (map<CParticle4D, int>::iterator it = cores.begin(); it != cores.end(); ++it)
	{
		numclusters = Max(numclusters, it->second + 1);
	}
	vector<vector<CParticle4D>> Q(numclusters);
	for (map<CParticle4D, int>::iterator it = cores.begin(); it != cores.end(); ++it)
	{
		Q[it->second].push_back(it->first);
		SetVoxel(S, it->first, it->second, dims);
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
				float dval = GetVoxel(D, p, 0.0f, dims);
				for (int n = 0; n < NumNeighborsHere; ++n)
				{
					CParticle4D q(p.m_X + XOffsetHere[n], p.m_Y + YOffsetHere[n]);
					if (GetVoxel(L, q, (unsigned char)0, dims) && GetVoxel(D, q, 0.0f, dims) < dval)
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
	//separableConvolution(D, D0, filter, CBE_ZERO, ndimL, dimsL);
	D = D0;

	map<CParticle4D,CParticle4D> coremap = SteepestGrowthHere(vfront, D, dimsL);
	
	vector<Vertex<CParticle4D>*> vertices = makeGraph(coremap, thres);
	map<CParticle4D,int> corelabel = clusterVertices(vertices, thres);
	//map<CParticle4D, int> surflabel = partitionSurface(vfront, corelabel, D, dimsL);

	vector<int> S(nvoxels, 0);
	partitionRegion(S, L, D, corelabel, minsize, ndimL, dimsL);

	if(nlhs >= 1)
	{
		plhs[0] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 2)
	{
		int dims[] = { corelabel.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (map<CParticle4D, int>::iterator it = corelabel.begin(); it != corelabel.end(); ++it, i++)
		{
			SetData2(F, i, 0, dims[0], dims[1], it->first.m_X);
			SetData2(F, i, 1, dims[0], dims[1], it->first.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], it->second);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		/*int dims[] = { factory.edges.size(), 5 };
		vector<float> F(dims[0] * dims[1]);
		for (int i = 0; i < dims[0]; ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], (float)factory.edges[i]->u->key.m_X);
			SetData2(F, i, 1, dims[0], dims[1], (float)factory.edges[i]->u->key.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], (float)factory.edges[i]->v->key.m_X);
			SetData2(F, i, 3, dims[0], dims[1], (float)factory.edges[i]->v->key.m_Y);
			SetData2(F, i, 4, dims[0], dims[1], (float)factory.edges[i]->w);
		}*/
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
		plhs[2] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(D, mxSINGLE_CLASS, ndimL, dimsL);
	}
	factory.Clean();
	mexUnlock();
}

