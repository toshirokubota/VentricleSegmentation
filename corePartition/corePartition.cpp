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

//Change these settings for another neighborhood system
const int NumNeighborsHere = 8; // 3 * NumNeighbors + 2;
int TOffsetArray[NumNeighborsHere];
int XOffsetArray[NumNeighborsHere];
int YOffsetArray[NumNeighborsHere];
int ZOffsetArray[NumNeighborsHere];

#define XOffsetHere XOffsetArray
#define YOffsetHere YOffsetArray
#define ZOffsetHere ZOffsetArray
#define TOffsetHere TOffsetArray
float OffsetLength[NumNeighborsHere];

void
initializeOffsetArrays()
{
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
set<CParticle4D>
steepestUpwardPropagation(CParticle4D p0, set<CParticle4D>& target, const vector<float>& D, const int* dims)
{
	set<CParticle4D> core;
	vector<CParticle4D> Q(1, p0);
	vector<unsigned char> M(D.size(), (unsigned char)0);
	SetVoxel(M, p0, (unsigned char)1, dims);

	while (Q.empty()==false)
	{
		vector<CParticle4D> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CParticle4D p = Q[i];
			if (target.find(p) != target.end())
			{
				core.insert(p);
				continue;
			}

			float dval1 = GetVoxel(D, p, 0.0f, dims);
			bool bmax = true;
			for (int j = 0; j<NumNeighborsHere; ++j)
			{
				CParticle4D q(p.m_X + XOffsetHere[j], p.m_Y + YOffsetHere[j]);
				if (GetVoxel(D, q, 0.0f, dims)>dval1)
				{
					vector<float> g = Gradient2(D, q, dims);
					float d = (g[0] * (q.m_X - p0.m_X) + g[1] * (q.m_Y - p0.m_Y));
					if (d > 0)
					{
						bmax = false;
						if (GetVoxel(M, q, (unsigned char)1, dims) == 0)
						{
							SetVoxel(M, q, (unsigned char)1, dims);
							Q2.push_back(q);
						}
					}
				}
			}
			if (bmax)
			{
				core.insert(p);
			}
		}
		Q = Q2;
	}
	return core;
}

/*
From each surface pixel (FRONT), find core pixels by going up straight in the distance map (D).
*/
set<CParticle4D>
locateMediaxAxis(vector<CParticle4D>& front,
		   const vector<float>& D,
		   const int* dims)
{
	set<CParticle4D> core;
	set<CParticle4D> target; //dummy
	for (vector<CParticle4D>::iterator it = front.begin(); it != front.end(); ++it)
	{
		set<CParticle4D> res = steepestUpwardPropagation(*it, target, D, dims);
		core.insert(res.begin(), res.end());
	}

	return core;
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
makeGraph(set<CParticle4D>& core, float thres)
{
	GraphFactory<CParticle4D>& factory = GraphFactory<CParticle4D>::GetInstance();
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
clusterVertices(vector<Vertex<CParticle4D>*>& vertices)
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

map<CParticle4D, int>
partitionSurface(vector<CParticle4D>& surface,
					map<CParticle4D, int>& label,
					set<CParticle4D>& coreset,
					float prop, //measure of dominance
					vector<float>& D, const int* dims)
{
	int numclusters = 0;
	for (map<CParticle4D, int>::iterator it = label.begin(); it != label.end(); ++it)
	{
		numclusters = Max(numclusters, it->second + 1);
	}
	map<CParticle4D, int> lmap;
	for (int i = 0; i < surface.size(); ++i)
	{
		CParticle4D p = surface[i];
		set<CParticle4D> res = steepestUpwardPropagation(p, coreset, D, dims);
		vector<int> count(numclusters, 0);
		for (set<CParticle4D>::iterator it = res.begin(); it != res.end(); ++it)
		{
			CParticle4D q = *it;
			count[label[q]]++;
		}
		vector<pair<int, int>> pairs(count.size()); 
		for (int j = 0; j < count.size(); ++j)
		{
			pairs[j] = pair<int, int>(count[j], j);
		}
		sort(pairs.begin(), pairs.end());
		if (pairs[pairs.size() - 1].first > pairs[pairs.size() - 2].first * prop)
		{
			lmap[p] = pairs[pairs.size() - 1].second;
		}
		else
		{
			lmap[p] = -1; //ambigous
		}
	}
	return lmap;
}

void
partitionRegion(vector<int>& S, //labeled image
				const vector<unsigned char>& L, //foreground
				map<CParticle4D, int>& cores, 
				map<CParticle4D,int>& surface,
				int minsize,
				int ndim, const int* dims)
{
	int numclusters = 0;
	for (map<CParticle4D, int>::iterator it = cores.begin(); it != cores.end(); ++it)
	{
		numclusters = Max(numclusters, it->second + 1);
		SetVoxel(S, it->first, it->second, dims);
	}
	vector<int> count(numclusters, 0);
	for (map<CParticle4D, int>::iterator it = surface.begin(); it != surface.end(); ++it)
	{
		count[it->second]++;
	}
	for (int i = 0; i < count.size(); ++i)
	{
		printf("count %d %d\n", i, count[i]);
	}
	vector<CParticle4D> Q;
	for (map<CParticle4D, int>::iterator it = cores.begin(); it != cores.end(); ++it)
	{
		if (count[it->second] >= minsize)
		{
			Q.push_back(it->first);
		}
	}
	for (int i = 0; i < dims[1]; ++i)
	{
		for (int j = 0; j < dims[0]; ++j)
		{
			CParticle4D p(j, i);
			if (GetVoxel(L, p, (unsigned char)0, dims))
			{
				float mind = std::numeric_limits<float>::infinity();
				CParticle4D q;
				for (int k = 0; k < Q.size(); ++k)
				{
					float d = Length2(p.m_X - Q[k].m_X, p.m_Y - Q[k].m_Y);
					if (d < mind)
					{
						mind = d;
						q = Q[k];
					}
				}
				if (mind < std::numeric_limits<float>::infinity())
				{
					SetVoxel(S, p, cores[q], dims);
				}
			}
		}
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
	float proportion = 0.5;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(proportion, prhs[2], classMode);
	}
	int minsize=1;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minsize, prhs[3], classMode);
	}
	printf("threshold = %f, proportion = %f, minsize=%d\n", thres, proportion, minsize);

	int nvoxels = numberOfElements(ndimL,dimsL);

	vector<CParticle4D> vfront = ExtactSurface(L, dimsL);
	printf("There are %d boundary voxels.\n", vfront.size());

	vector<float> D(nvoxels, 0);
	vector<unsigned char> iL(nvoxels, 0);
	for (int i = 0; i<nvoxels; ++i)
	{
		iL[i] = L[i] ? 0 : 1;
	}
	vector<float> v(2, 1.0f);
	DistanceTransformEuclidF(D, iL, v, ndimL, dimsL);
	set<CParticle4D> coreset = locateMediaxAxis(vfront, D, dimsL);

	GraphFactory<CParticle4D>& factory = GraphFactory<CParticle4D>::GetInstance();
	vector<Vertex<CParticle4D>*> vertices = makeGraph(coreset, thres);
	map<CParticle4D, int> corelabel = clusterVertices(vertices);
	map<CParticle4D, int> surflabel = partitionSurface(vfront, corelabel, coreset, proportion, D, dimsL);
	vector<int> S(nvoxels, 0);
	partitionRegion(S, L, surflabel, surflabel, minsize, ndimL, dimsL);

	if (nlhs >= 1)
	{
		int dims[] = { vertices.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < vertices.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], vertices[i]->key.m_X);
			SetData2(F, i, 1, dims[0], dims[1], vertices[i]->key.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], corelabel[vertices[i]->key]);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { surflabel.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (map<CParticle4D,int>::iterator it = surflabel.begin(); it != surflabel.end(); ++it, i++)
		{
			SetData2(F, i, 0, dims[0], dims[1], it->first.m_X);
			SetData2(F, i, 1, dims[0], dims[1], it->first.m_Y);
			SetData2(F, i, 2, dims[0], dims[1], it->second);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(D, mxSINGLE_CLASS, ndimL, dimsL);
	}
	mexUnlock();
}

