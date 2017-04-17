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
#include <szMiscOperations.h>
#include <DisjointSet.h>
#include <szParticle4D.h>
#include <szConnectedComponent.h>
#include <CoreParticle.h>
#include <NeighborhoodND.h>
#include <CoreParticleUtil.h>
#include <CoreParticlePropagation.h>
#include <CoreParticleUtilTemplate.h>

float
dotProduct(float x1, float y1, float z1, float t1, float x2, float y2, float z2, float t2)
{
	return x1*x2 + y1*y2 + z1*z2 + t1*t2;
}

float
length(float x, float y, float z, float t)
{
	return sqrt(x*x + y*y + z*z + t*t);
}

CoreParticle* trace2Medial(CoreParticle* p, CoreParticle* q, int ndim)
{
	float dir[4] = { q->x - p->x, q->y - p->y, q->z - p->z, q->t - p->t };
	while (medialParticle(q, ndim) == false)
	{
		float maxdp = -std::numeric_limits<float>::infinity();
		CoreParticle* chosen = NULL;
		for (set<CoreParticle*>::iterator it = q->descendents.begin(); it != q->descendents.end(); ++it)
		{
			CoreParticle* r = *it;
			float dp = dotProduct(r->x - q->x, r->y - q->y, r->z - q->z, r->t - q->t, dir[0], dir[1], dir[2], dir[3]);
			float len = length(r->x - q->x, r->y - q->y, r->z - q->z, r->t - q->t);
			dp /= len;
			if (dp > maxdp)
			{
				maxdp = dp;
				chosen = r;
			}
		}
		p = q;
		q = chosen;
	}
	return q;
}

vector<Vertex<CoreParticle*>*>
pruneEdges(vector<Edge<CoreParticle*>*>& edges) {
	set<Vertex<CoreParticle*>*> sv;
	set<Edge<CoreParticle*>*> se;
	for (int i = 0; i < edges.size(); ++i)
	{
		sv.insert(edges[i]->u);
		sv.insert(edges[i]->v);
		se.insert(edges[i]);
		//add the other direction, too.
		Edge<CoreParticle*>* e2 = edges[i]->v->findEdge(edges[i]->u);
		if (e2 != NULL)
		{
			se.insert(e2);
		}
	}
	se.insert(edges.begin(), edges.end());
	for (set<Vertex<CoreParticle*>*>::iterator it = sv.begin(); it != sv.end(); ++it)
	{
		Vertex<CoreParticle*>* v = *it;
		for (vector<Edge<CoreParticle*>*>::iterator jt = v->aList.end() - 1; jt >= v->aList.begin(); jt--)
		{
			if (find(se.begin(), se.end(), *jt) == se.end())
			{
				v->aList.erase(jt);
			}
		}
	}
	vector<Vertex<CoreParticle*>*> vv;
	vv.insert(vv.begin(), sv.begin(), sv.end());

	return vv;
}

bool
DFS(Vertex<CoreParticle*>* u, set<Vertex<CoreParticle*>*>& cover)
{
	u->color = Black;
	bool bFound = cover.find(u) != cover.end();
	for (int i = 0; i < u->aList.size(); ++i)
	{
		Vertex<CoreParticle*>* v = u->aList[i]->v;
		if (v->color == White)
		{
			if (DFS(v, cover))
			{
				bFound = true;
			}
		}
	}
	if (bFound)
	{
		cover.insert(u);
	}

	return bFound;
}

/*
Traverse the tree from one of vertex in the cover with DFS.
If it reaches another vertex in the cover, add all the intermediates to the cover.
*/
set<Vertex<CoreParticle*>*> 
collectInBetween(set<Vertex<CoreParticle*>*>& s, vector<Vertex<CoreParticle*>*>& vertices)
{
	set<Vertex<CoreParticle*>*> cover;
	cover.insert(s.begin(), s.end());

	for (vector<Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it)
	{
		(*it)->color = White;
	}
	Vertex<CoreParticle*>* o = *cover.begin();
	DFS(o, cover);

	return cover;
}

CoreParticle* findTheOtherSide(CoreParticle* p, CoreParticle* q,
	vector<CoreParticle*>& mp, int ndim, const int* dims,
	bool bEightNeighbor=false)
{
	float dir[4] = { q->x - p->x, q->y - p->y, q->z - p->z, q->t - p->t };
	while (surfaceParticle(q, ndim, bEightNeighbor) == false)
	{

		float maxdp = 0;
		CoreParticle* chosen = NULL;
		for (vector<CoreParticle*>::iterator it = q->neighbors.begin(); it != q->neighbors.end(); ++it)
		{
			CoreParticle* r = *it;
			float dp = dotProduct(r->x - q->x, r->y - q->y, r->z - q->z, r->t - q->t, dir[0], dir[1], dir[2], dir[3]);
			float len = length(r->x - q->x, r->y - q->y, r->z - q->z, r->t - q->t);
			dp /= len;
			if (dp > maxdp)
			{
				maxdp = dp;
				chosen = r;
			}
		}
		q = chosen;
	}
	return q;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [S] = particleGrowth(L, [thres mincore])");
		return;
	}

	//load figure-ground segmentation
	int ndimL;
	const int* dimsL;
	mxClassID classL;
	vector<unsigned char> L;
	LoadData(L, prhs[0], classL, ndimL, &dimsL);

	//int nclusters = 2;
	int cutoff = 3;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(cutoff, prhs[1], classMode);
	}
	bool _EightNeighbor = false;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		int value;
		ReadScalar(value, prhs[2], classMode);
		_EightNeighbor = value > 0 ? true : false;
	}

	int nvoxels = numberOfElements(ndimL, dimsL);

	vector<CoreParticle*> mp = generateParticleMap(L, ndimL, dimsL);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndimL, dimsL);
	vector<int> S(nvoxels, 0);
	propagateParticles(particles, mp, S, ndimL, dimsL);
	vector<int> S2(nvoxels, 0);
	vector<Edge<CoreParticle*>*> mst = makeGraphStructure(mp, S2, ndimL, dimsL);
	vector<Vertex<CoreParticle*>*> vertices = pruneEdges(mst);
	map<CoreParticle*, Vertex<CoreParticle*>*> vmap;
	for (int i = 0; i < vertices.size(); ++i)
	{
		vmap[vertices[i]->key] = vertices[i];
	}

	vector<CoreParticle*> surface;
	vector<CoreParticle*> inflection;
	vector<CoreParticle*> medial;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (surfaceParticle(particles[i], ndimL, false))
		{
			surface.push_back(particles[i]);
		}
		if (inflectionParticle(particles[i], ndimL))
		{
			inflection.push_back(particles[i]);
		}
		if (medialParticle(particles[i], ndimL))
		{
			medial.push_back(particles[i]);
		}
	}

	//find the best cutting plane from each inflection point and extract a part separated from the main body.
	for (int i = 0; i < inflection.size(); ++i)
	{
		CoreParticle* p = inflection[i];
		set<Vertex<CoreParticle*>*> medaix;
		for (set<CoreParticle*>::iterator it = p->descendents.begin(); it != p->descendents.end(); ++it)
		{
			CoreParticle* z = trace2Medial(p, *it, ndimL);
			if (find(vmap.begin(), vmap.end(), z) == vmap.end())
			{
				printf("Warning: medial axis particle not in MST");
			}
			else
			{
				medaix.insert(vmap[z]);
			}
		}
		if (medaix.size() > 1)
		{
			printf("Warning: there are only %d medial axis intersections.", medaix.size());
		}

		set<Vertex<CoreParticle*>*> inbetweens = collectInBetween(medaix, vertices);

		for (set<Vertex<CoreParticle*>*>::iterator it = inbetweens.begin(); it != inbetweens.end(); ++it)
		{
			CoreParticle* q = findTheOtherSide(p, (*it)->key, mp, ndimL, dimsL, _EightNeighbor);
			vector<CoreParticle*> boundary = cutBoundary(p, q, surface);
		}
	}

	if (nlhs >= 1)
	{
		int dims[] = { particles.size(), 10 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->z);
			SetData2(F, i, 3, dims[0], dims[1], p->t);
			SetData2(F, i, 4, dims[0], dims[1], p->label);
			SetData2(F, i, 5, dims[0], dims[1], p->gen);
			SetData2(F, i, 6, dims[0], dims[1], p->id);
			SetData2(F, i, 7, dims[0], dims[1], surfaceParticle(p, ndimL) ? 1 : 0);
			SetData2(F, i, 8, dims[0], dims[1], medialParticle(p, ndimL) ? 1 : 0);
			SetData2(F, i, 9, dims[0], dims[1], inflectionParticle(p, ndimL) ? (int)p->descendents.size() : 0);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { adj.size(), 2 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (set<pair<CoreParticle*, CoreParticle*>>::iterator it = adj.begin(); it != adj.end(); ++it)
		{
			pair<CoreParticle*, CoreParticle*> pr = *it;
			SetData2(F, i, 0, dims[0], dims[1], pr.first->id);
			SetData2(F, i, 1, dims[0], dims[1], pr.second->id);
			i++;
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(S2, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 5)
	{
		vector<pair<int, int>> pairs;
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			for (set<CoreParticle*>::iterator it = p->descendents.begin(); it != p->descendents.end(); ++it)
			{
				pairs.push_back(pair<int, int>(p->id, (*it)->id));
			}
		}
		int dims[] = { pairs.size(), 2 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < pairs.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], pairs[i].first);
			SetData2(F, i, 1, dims[0], dims[1], pairs[i].second);
		}
		plhs[4] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	GraphFactory<CoreParticle*>::GetInstance().Clean();
	NeighborhoodFactory::getInstance().clean();
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}

