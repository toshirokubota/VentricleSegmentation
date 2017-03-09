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
#include <szDistanceTransform.h>
#include <szConnectedComponent.h>
#include <szLocalExtrema.h>
#include <szDistanceTransformNonIsotropic.h>
#include <Graph.h>
#include <GraphFactory.h>

bool _EightNeighbor = false;

struct CoreParticle
{
	CoreParticle(int x0 = 0, int y0 = 0, int z0 = 0, int t0 = 0, int gen0 = 0, int lb = 0, int id0 = 0)
	{
		x = x0;
		y = y0;
		gen = gen0;
		label = lb;
		id = id0;
		z = z0;
		t = t0;
		value = 0;
		pi = NULL;
		selected = false;
	}
	int x;
	int y;
	int z;
	int t;
	int label;
	int gen; //generation
	int id;
	float value; //generic value
	bool selected;
	set<CoreParticle*> ascendents;
	set<CoreParticle*> descendents;
	vector<CoreParticle*> neighbors;
	CoreParticle* pi; //path from the core
};

struct CoreParticleFactory
{
public:
	static CoreParticleFactory& getInstance()
	{
		static CoreParticleFactory instance;
		return instance;
	}
	CoreParticle* makeParticle(int x = 0, int y = 0, int z = 0, int t = 0, int gen = 0, int lb = 0)
	{
		CoreParticle* particle = new CoreParticle(x, y, z, t, gen, lb, _id++);
		particles.push_back(particle);
		return particle;
	}
	void clean()
	{
		for (int i = 0; i < particles.size(); ++i)
		{
			delete particles[i];
		}
		particles.clear();
		_id = 0;
	}
	vector<CoreParticle*> particles;
private:
	int _id;
	CoreParticleFactory()
	{
		_id = 0;
	}
	~CoreParticleFactory()
	{
		clean();
	}
	CoreParticleFactory(CoreParticleFactory& f){}
	CoreParticleFactory operator=(CoreParticleFactory& f){}
};

struct NeighborhoodFactory
{
public:
	static NeighborhoodFactory& getInstance(int n = 0)
	{
		static NeighborhoodFactory instance(n);
		if (n > 0 && n != instance.ndim)
		{
			instance.neighbor4 = MakeFourNeighborhood(n);
			instance.neighbor8 = MakeEightNeighborhood(n);
			instance.ndim = n;
		}
		return instance;
	}
	vector<vector<int>> neighbor4;
	vector<vector<int>> neighbor8;
private:
	int ndim;
	NeighborhoodFactory(int n)
	{
		ndim = n;
		if (ndim > 0)
		{
			printf("I am here.");

			neighbor4 = MakeFourNeighborhood(n);
			neighbor8 = MakeEightNeighborhood(n);
		}
	}
	~NeighborhoodFactory()
	{
	}
	NeighborhoodFactory(NeighborhoodFactory& f){}
	NeighborhoodFactory operator=(NeighborhoodFactory& f){}
};


CoreParticle*
coreParticleNdim(vector<int>& loc, int ndim)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	switch (ndim)
	{
	case 1:
		return factory.makeParticle(loc[0]);
	case 2:
		return factory.makeParticle(loc[0], loc[1]);
	case 3:
		return factory.makeParticle(loc[0], loc[1], loc[2]);
	case 4:
		return factory.makeParticle(loc[0], loc[1], loc[2], loc[3]);
	default:
		mexErrMsgTxt("particleNdim: unsupported number of dimensions. It has to be between 1 and 4.");
		return NULL;
	}
}

vector<int>
coreParticle2Index(CoreParticle*  p, int ndim)
{
	vector<int> idx(ndim);
	if (ndim == 1)
	{
		idx[0] = p->x;
	}
	else if (ndim == 2)
	{
		idx[0] = p->x;
		idx[1] = p->y;
	}
	else if (ndim == 3)
	{
		idx[0] = p->x;
		idx[1] = p->y;
		idx[2] = p->z;
	}
	else if (ndim == 4)
	{
		idx[0] = p->x;
		idx[1] = p->y;
		idx[2] = p->z;
		idx[3] = p->t;
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
	}
	return idx;
}

template<class T>
bool
SetVoxel(vector<T>& A,
const CoreParticle* p,
T value,
int ndim,
const int* dims)
{
	if (ndim == 1)
	{
		return SetData(A, p->x, value);
	}
	else if (ndim == 2)
	{
		return SetData2(A, p->x, p->y, dims[0], dims[1], value);
	}
	else if (ndim == 3)
	{
		return SetData3(A, p->x, p->y, p->z, dims[0], dims[1], dims[2], value);
	}
	else if (ndim == 4)
	{
		return SetData4(A, p->x, p->y, p->z, p->t, dims[0], dims[1], dims[2], dims[3], value);
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
		return false;
	}
}

template<class T>
T
GetVoxel(const vector<T>& A,
const CoreParticle* p,
T defaultValue,
int ndim,
const int* dims)
{
	if (ndim == 1)
	{
		return GetData(A, p->x, defaultValue);
	}
	else if (ndim == 2)
	{
		return GetData2(A, p->x, p->y, dims[0], dims[1], defaultValue);
	}
	else if (ndim == 3)
	{
		return GetData3(A, p->x, p->y, p->z, dims[0], dims[1], dims[2], defaultValue);
	}
	else if (ndim == 4)
	{
		return GetData4(A, p->x, p->y, p->z, p->t, dims[0], dims[1], dims[2], dims[3], defaultValue);
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
		return defaultValue;
	}
}


vector<CoreParticle*>
generateParticleMap(vector<unsigned char>& L,
int ndim,
const int* dims)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	for (int i = 0; i < L.size(); ++i)
	{
		if (L[i])
		{
			vector<int> sub = Ind2Sub(i, ndim, dims);
			CoreParticle* p = coreParticleNdim(sub, ndim);
			mp[i] = p;
		}
	}
	return mp;
}

vector<CoreParticle*>
setupParticleNeighbors(vector<CoreParticle*>& mp,
int ndim,
const int* dims)
{
	vector<CoreParticle*> particles;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i] != NULL)
		{
			particles.push_back(mp[i]);
		}
	}
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4;
	if (_EightNeighbor)
	{
		nbh = NeighborhoodFactory::getInstance(ndim).neighbor8; //MakeEightNeighborhood(ndim, dims);
	}

	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		vector<int> sub = coreParticle2Index(p, ndim);
		for (int n = 0; n < nbh.size(); ++n)
		{
			if (NeighborCheck(sub.begin(), nbh[n].begin(), ndim, dims))
			{
				int idx = Sub2Ind(sub, nbh[n], ndim, dims);
				if (mp[idx] != NULL)
				{
					p->neighbors.push_back(mp[idx]);
				}
			}
		}
	}
	return particles;
}

/*
Cluster particles using disjoint set.
*/
vector<vector<CoreParticle*>>
clusterParticles(vector<CoreParticle*>& particles)
{
	set<CoreParticle*> S;
	vector<Node<CoreParticle*>*> nodes;
	map<CoreParticle*, int> imap;
	for (int i = 0; i < particles.size(); ++i)
	{
		Node<CoreParticle*>* n = makeset(particles[i]);
		nodes.push_back(n);
		imap[particles[i]] = i;
		S.insert(particles[i]);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		for (int j = 0; j < p->neighbors.size(); ++j)
		{
			CoreParticle* q = p->neighbors[j];
			if (S.find(q) != S.end())
			{
				merge(nodes[i], nodes[imap[q]]);
			}
		}
	}
	vector<Node<CoreParticle*>*> rep = clusters(nodes);
	vector<vector<CoreParticle*>> group(rep.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = distance(rep.begin(), find(rep.begin(), rep.end(), findset(nodes[i])));
		group[k].push_back(nodes[i]->key);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return group;
}

void
propagateDescendency(vector<CoreParticle*>& P)
{
	vector<CoreParticle*> Q;
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->descendents.size()>0)
		{
			Q.push_back(P[i]);
		}
	}
	while (Q.empty() == false)
	{
		set<CoreParticle*> S;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			for (int j = 0; j < p->neighbors.size(); ++j)
			{
				CoreParticle* q = p->neighbors[j];
				if (p->gen == q->gen)
				{
					if (q->descendents.size() == 0)
					{
						q->descendents.insert(p);
						p->ascendents.insert(q);
						S.insert(q);
					}
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
}

/*
from each surface particle, move towrad the center and establish ascendent/descendent relation.
*/
vector<CoreParticle*>
propagateParticles(
vector<CoreParticle*>& particles,
vector<CoreParticle*>& mp,
vector<int>& S,
int ndim,
const int* dims)
{
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4;
	if (_EightNeighbor)
	{
		nbh = NeighborhoodFactory::getInstance(ndim).neighbor8; //MakeEightNeighborhood(ndim, dims);
	}

	vector<CoreParticle*> Q;
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		if (p->neighbors.size() < nbh.size())
		{
			Q.push_back(p);
		}
	}
	int gen = 1;
	while (Q.empty() == false)
	{
		//vector<CoreParticle*> core;
		for (int i = 0; i < Q.size(); ++i)
		{
			SetVoxel(S, Q[i], gen, ndim, dims);
			Q[i]->gen = gen;
		}

		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			CoreParticle* p2 = NULL;
			for (int n = 0; n < p->neighbors.size(); ++n)
			{
				CoreParticle* q = p->neighbors[n];
				int sval = GetVoxel(S, q, 0, ndim, dims);
				if (sval == 0)
				{
					Q2.insert(q);
					p->descendents.insert(q);
					q->ascendents.insert(p);
				}
			}
		}
		propagateDescendency(Q);

		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
		gen++;
	}
	return particles;
}

vector<CoreParticle*>
localGenMaximum(vector<CoreParticle*>& pmap, int ndim, const int* dims)
{
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeEightNeighborhood(ndim, dims);
	vector<CoreParticle*> trace;

	for (int i = 0; i < pmap.size(); ++i)
	{
		CoreParticle* p = pmap[i];
		if (p == NULL) continue;

		bool bMax = true;
		for (int n = 0; n < p->neighbors.size(); ++n)
		{
			CoreParticle* q = p->neighbors[n];
			if (p->gen < q->gen)
			{
				bMax = false;
			}
		}
		if (bMax)
		{
			p->selected = true;
			trace.push_back(p);
		}
	}
	return trace;
}


vector<CoreParticle*>
removeSpurious(vector<CoreParticle*>& cores, float thres)
{
	vector<vector<CoreParticle*>> cl = clusterParticles(cores);
	vector<CoreParticle*> result;
	vector<CoreParticle*> removed;
	for (int i = 0; i < cl.size(); ++i)
	{
		set<CoreParticle*> up;
		set<CoreParticle*> down;
		for (int j = 0; j < cl[i].size(); ++j)
		{
			CoreParticle* p = cl[i][j];
			for (int k = 0; k < p->neighbors.size(); ++k)
			{
				CoreParticle* q = p->neighbors[k];
				if (q->selected == false)
				{
					if (q->gen >= p->gen)
					{
						up.insert(q);
					}
					else
					{
						down.insert(q);
					}
				}
			}
		}
		if ((float)down.size() / ((float)up.size() + down.size()) > thres)
		{
			result.insert(result.end(), cl[i].begin(), cl[i].end());
		}
		else
		{
			removed.insert(removed.end(), cl[i].begin(), cl[i].end());
		}
	}
	for (int i = 0; i < removed.size(); ++i)
	{
		cores[i]->selected = false;
		CoreParticle* p = removed[i];
		for (int j = 0; j < p->neighbors.size(); ++j)
		{
			CoreParticle* q = p->neighbors[j];
			if (p->gen == q->gen && q->descendents.empty() == false)
			{
				//p->descendents.insert(q);
				//q->ascendents.insert(p);
			}
		}
	}

	return result;
}

void
colorParticles(vector<CoreParticle*>& cores, int thres, int ndim, const int* dims)
{
	vector<vector<CoreParticle*>> cl = clusterParticles(cores);
	vector<CoreParticle*> core;
	//int maxgen = 0;
	for (int i = 0; i < cl.size(); ++i)
	{
		if (cl[i][0]->gen <= thres)
		{
			core.insert(core.begin(), cl[i].begin(), cl[i].end());
			//maxgen = cl[i][0]->gen;
		}
	}
	vector<CoreParticle*> Q = core;
	while (Q.empty() == false)
	{
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
			//for (vector<CoreParticle*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
			{
				CoreParticle* q = *it;
				//if (q->gen < p->gen)
				{
					q->label = 1;
					Q2.insert(q);
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
	}
}

void
colorParticles2(vector<CoreParticle*>& cores, int niter, int ndim, const int* dims)
{
	vector<vector<CoreParticle*>> cl = clusterParticles(cores);
	vector<CoreParticle*> core;
	int maxgen = 0;
	for (int i = 0; i < cl.size(); ++i)
	{
		if (cl[i][0]->gen > maxgen)
		{
			core = cl[i];
			maxgen = cl[i][0]->gen;
		}
	}
	vector<CoreParticle*> Q = core;
	for (int iter = 0; iter < niter; ++iter)
	{
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			//for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
			for (vector<CoreParticle*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
			{
				CoreParticle* q = *it;
				//if (q->gen < p->gen)
				{
					q->label = 1;
					Q2.insert(q);
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
	}
}

vector<vector<CoreParticle*>>
mergeClusters(vector<CoreParticle*>& core, vector<CoreParticle*>& mp, vector<int>& T, int ndim, const int* dims)
{
	vector<vector<CoreParticle*>> groups = clusterParticles(core);
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			mp[i]->value = std::numeric_limits<float>::infinity();
			mp[i]->label = -1;
		}
	}
	for (int i = 0; i < groups.size(); ++i)
	{
		for (int j = 0; j < groups[i].size(); ++j)
		{
			groups[i][j]->label = i;
			groups[i][j]->value = 0.0f;
		}
	}
	vector<CoreParticle*> Q = core;
	while (Q.empty() == false)
	{
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			SetVoxel(T, p, p->label + 1, ndim, dims);
			for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
			{
				CoreParticle* q = *it;
				if (q->label < 0)
				{
					q->label = p->label;
					q->value = p->value + 1;
					Q2.insert(q);
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
	}
	//find adjacent regions
	vector<vector<bool>> adj;
	vector<vector<float>> dist;
	for (int i = 0; i < groups.size(); ++i)
	{
		adj.push_back(vector<bool>(groups.size(), false));
		dist.push_back(vector<float>(groups.size(), std::numeric_limits<float>::infinity()));
	}
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			CoreParticle* p = mp[i];
			for (int j = 0; j < p->neighbors.size(); ++j)
			{
				CoreParticle* q = p->neighbors[j];
				if (p->label != q->label)
				{
					adj[p->label][q->label] = true;
					adj[q->label][p->label] = true;
					dist[q->label][p->label] = Min(q->value, dist[q->label][p->label]);
					dist[p->label][q->label] = Min(p->value, dist[p->label][q->label]);
				}
			}
		}
	}

	//relable based on the bonding measure
	vector<Node<int>*> nodes;
	for (int i = 0; i < groups.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for (int i = 0; i < groups.size(); ++i)
	{
		for (int j = i + 1; j < groups.size(); ++j)
		{
			if (adj[i][j])
			{
				int dval = (int)(dist[j][i] + dist[i][j]);
				printf("%d vs. %d: gen=%d, %d, dist=%d, %d\n",
					i, j, groups[i][0]->gen, groups[j][0]->gen, (int)dist[i][j], (int)dist[j][i]);
				if (groups[i][0]->gen > dval && groups[j][0]->gen > dval)
				{
					merge(nodes[i], nodes[j]);
				}
			}
		}
	}
	vector<Node<int>*> reps = clusters(nodes);
	vector<vector<CoreParticle*>> result(reps.size());
	for (int i = 0; i < groups.size(); ++i)
	{
		int k = distance(reps.begin(), find(reps.begin(), reps.end(), findset(nodes[i])));
		for (int j = 0; j < groups[i].size(); ++j)
		{
			result[k].insert(result[k].end(), groups[i].begin(), groups[i].end());
		}
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return result;
}

void
colorParticles3(vector<CoreParticle*>& cores, vector<CoreParticle*>& mp, vector<int>& T, int ndim, const int* dims)
{
	vector<vector<CoreParticle*>> cl = mergeClusters(cores, mp, T, ndim, dims);
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			mp[i]->label = 0;
		}
	}
	set<int> sgen;
	for (int i = 0; i < cl.size(); ++i)
	{
		sgen.insert(cl[i][0]->gen);
		for (int j = 0; j < cl[i].size(); ++j)
		{
			cl[i][j]->label = i+1;
		}
	}
	vector<int> vgen;
	vgen.insert(vgen.begin(), sgen.begin(), sgen.end());
	sort(vgen.begin(), vgen.end());
	for (int ig = 0; ig < vgen.size(); ++ig)
	{
		int gval = vgen[ig];
		vector<CoreParticle*> Q;
		for (int i = 0; i < cl.size(); ++i)
		{
			if (cl[i][0]->gen == gval)
			{
				Q.insert(Q.end(), cl[i].begin(), cl[i].end());
			}
		}
		while (Q.empty()==false)
		{
			set<CoreParticle*> Q2;
			for (int i = 0; i < Q.size(); ++i)
			{
				CoreParticle* p = Q[i];
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				//for (vector<CoreParticle*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->label <= 0)
					{
						q->label = p->label;
						Q2.insert(q);
					}
				}
			}
			Q.clear();
			Q.insert(Q.end(), Q2.begin(), Q2.end());
		}
	}
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
	float cutoff = 0.5f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(cutoff, prhs[1], classMode);
	}
	//int niter = 5;
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

	vector<CoreParticle*> trace = localGenMaximum(mp, ndimL, dimsL);
	trace = removeSpurious(trace, cutoff);

	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
	}
	vector<int> T2(mp.size(), 0);
	colorParticles3(trace, mp, T2, ndimL, dimsL);
	vector<int> T(mp.size(), 0);
	for (int i = 0; i <	mp.size(); ++i)
	{
		if (mp[i] && mp[i]->label)
		{
			T[i] = mp[i]->label;
		}
	}

	if (nlhs >= 1)
	{
		int dims[] = { particles.size(), 7 };
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
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { trace.size(), 7 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < trace.size(); ++i)
		{
			CoreParticle* p = trace[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->z);
			SetData2(F, i, 3, dims[0], dims[1], p->t);
			SetData2(F, i, 4, dims[0], dims[1], p->label);
			SetData2(F, i, 5, dims[0], dims[1], p->gen);
			SetData2(F, i, 6, dims[0], dims[1], p->id);
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(T2, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(T, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 5)
	{
		plhs[4] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}

