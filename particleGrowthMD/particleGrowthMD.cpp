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

struct CoreParticle 
{
	CoreParticle(int x0 = 0, int y0 = 0, int z0=0, int t0=0, int gen0 = 0, int lb = 0, int id0 = 0)
	{
		x = x0;
		y = y0;
		gen = gen0;
		label = lb;
		id = id0;
		z = z0;
		t = t0;
		//corelink = NULL;
	}
	int x;
	int y;
	int z;
	int t;
	int label;
	int gen; //generation
	int id;
	set<CoreParticle*> ascendents;
	set<CoreParticle*> descendents;
	map<CoreParticle*, int> trace_time;
	//CoreParticle* corelink; //used to link cores
};

struct CoreParticleFactory
{
public:
	static CoreParticleFactory& getInstance()
	{
		static CoreParticleFactory instance;
		return instance;
	}
	CoreParticle* makeParticle(int x=0, int y=0, int z=0, int t=0, int gen=0, int lb = 0)
	{
		CoreParticle* particle = new CoreParticle(x, y, z, t, gen, lb, _id++);
		particles.push_back(particle);
		return particle;
	}
	void clean()
	{
		for (int i = 0; i<particles.size(); ++i)
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

struct cluster_pair;
struct core_cluster
{
	core_cluster(int id0)
	{
		id = id0;
		width = 0;
	}
	CParticle4D centroid() {
		CParticle4D c;
		for (int i = 0; i < particles.size(); ++i)
		{
			c.m_X += particles[i]->x;
			c.m_Y += particles[i]->y;
			c.m_Z += particles[i]->z;
			c.m_T += particles[i]->t;
		}
		c.m_X /= particles.size();
		c.m_Y /= particles.size();
		c.m_Z /= particles.size();
		c.m_T /= particles.size();
		return c;
	}
	vector<CoreParticle*> particles;
	vector<cluster_pair*> pairs;
	int size() { return particles.size(); }
	int width;
	int id;
};

struct cluster_pair
{
	cluster_pair() {
		u = NULL;
		v = NULL;
		cut = std::numeric_limits<float>::infinity();
		length = 0;
	}
	core_cluster* u;
	core_cluster* v;
	float cut;
	float length;
	float weight() {
		return (float)length * sqrt(u->width*v->width) / (float)cut * Min(u->size(), v->size()); // / Max(u->size(), v->size());
	}
};


struct CoreClusterFactory
{
public:
	static CoreClusterFactory& getInstance()
	{
		static CoreClusterFactory instance;
		return instance;
	}
	core_cluster* makeCoreCluster()
	{
		core_cluster* c = new core_cluster(_id++);
		clusters.push_back(c);
		return c;
	}
	cluster_pair* makeClusterPair(core_cluster* a, core_cluster* b)
	{
		cluster_pair* p = new cluster_pair();
		p->u = a;
		p->v = b;
		pairs.push_back(p);
		return p;
	}
	void clean()
	{
		for (int i = 0; i<clusters.size(); ++i)
		{
			delete clusters[i];
		}
		clusters.clear();
		for (int i = 0; i<pairs.size(); ++i)
		{
			delete pairs[i];
		}
		pairs.clear();
		_id = 0;
	}
	vector<core_cluster*> clusters;
	vector<cluster_pair*> pairs;
private:
	int _id;
	CoreClusterFactory()
	{
		_id = 0;
	}
	~CoreClusterFactory()
	{
		clean();
	}
	CoreClusterFactory(CoreClusterFactory& f){}
	CoreClusterFactory operator=(CoreClusterFactory& f){}
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

core_cluster*
merge_cores(core_cluster* a, core_cluster* b)
{
	CoreClusterFactory& factory = CoreClusterFactory::getInstance();
	core_cluster* c = factory.makeCoreCluster();
	c->particles = a->particles;
	c->particles.insert(c->particles.end(), b->particles.begin(), b->particles.end());
	c->width = Max(a->width, b->width);  //Is this right??
	for (int i = 0; i < a->pairs.size(); ++i)
	{
		if (a->pairs[i]->v != b)
		{
			a->pairs[i]->u = c;
			core_cluster* v = a->pairs[i]->v;
			for (int j = 0; j < v->pairs.size(); ++j)
			{
				if (v->pairs[j]->v == a)
				{
					v->pairs[j]->v = c;
				}
			}
			c->pairs.push_back(a->pairs[i]);
		}
	}
	for (int i = 0; i < b->pairs.size(); ++i)
	{
		if (b->pairs[i]->v != a)
		{
			b->pairs[i]->u = c;
			core_cluster* v = b->pairs[i]->v;
			for (int j = 0; j < v->pairs.size(); ++j)
			{
				if (v->pairs[j]->v == b)
				{
					v->pairs[j]->v = c;
				}
			}
			c->pairs.push_back(b->pairs[i]);
		}
	}
	return c;
}

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
const CParticle4D& p,
T value,
int ndim,
const int* dims)
{
	if (ndim == 1)
	{
		return SetData(A, p.m_X, value);
	}
	else if (ndim == 2)
	{
		return SetData2(A, p.m_X, p.m_Y, dims[0], dims[1], value);
	}
	else if (ndim == 3)
	{
		return SetData3(A, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], value);
	}
	else if (ndim == 4)
	{
		return SetData4(A, p.m_X, p.m_Y, p.m_Z, p.m_T, dims[0], dims[1], dims[2], dims[3], value);
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
const CParticle4D& p,
T defaultValue,
int ndim,
const int* dims)
{
	if (ndim == 1)
	{
		return GetData(A, p.m_X, defaultValue);
	}
	else if (ndim == 2)
	{
		return GetData2(A, p.m_X, p.m_Y, dims[0], dims[1], defaultValue);
	}
	else if (ndim == 3)
	{
		return GetData3(A, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], defaultValue);
	}
	else if (ndim == 4)
	{
		return GetData4(A, p.m_X, p.m_Y, p.m_Z, p.m_T, dims[0], dims[1], dims[2], dims[3], defaultValue);
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
		return defaultValue;
	}
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

/*
onSurface3 with a general neighborhood
*/
template<class T>
bool
onSurface(const vector<T>& A, const vector<int>& sub, const vector<vector<int>>& nbh, int ndim, const int* dims)
{
	for (int m = 0; m < nbh.size(); ++m)
	{
		if (NeighborCheck(sub.begin(), nbh[m].begin(), ndim, dims) && A[Sub2Ind(sub, nbh[m], ndim, dims)] == 0)
		{
			return true;
		}
	}
	return false;
}

vector<CoreParticle*>
ExtactSurface(const vector<unsigned char>& L, 
			  int ndim,
			  const int* dims)
{
	int nvoxels = numberOfElements(ndim, dims);
	vector<CoreParticle*> front;
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeEightNeighborhood(ndim, dims);
	for (int i = 0; i < nvoxels; ++i)
	{
		if (L[i])
		{
			vector<int> sub = Ind2Sub(i, ndim, dims);
			if (onSurface(L, sub, nbh, ndim, dims))
			{
				front.push_back(coreParticleNdim(sub, ndim));
			}
		}
	}
	return front;
}

float 
Length(float x, float y, float z=0, float t=0)
{
	return sqrt(x*x + y*y + z*z + t*t);
}

float
Length(vector<int>& idx, int ndim)
{
	float sum = 0;
	for (int i = 0; i < ndim; ++i)
	{
		sum += idx[i] * idx[i];
	}
	return sqrt(sum);
}

float
Distance(CoreParticle* p, CoreParticle* q)
{
	return Length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
}

vector<CoreParticle*>
pickCore(vector<CoreParticle*>& particles)
{
	vector<CoreParticle*> core;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i]->descendents.empty())
		{
			core.push_back(particles[i]);
		}
	}
	return core;
}

/*
Cluster particles using disjoint set.
*/
vector<vector<CoreParticle*>>
clusterParticles(vector<CoreParticle*>& particles, float thres)
{
	vector<Node<CoreParticle*>*> nodes;
	for (int i = 0; i < particles.size(); ++i)
	{
		nodes.push_back(makeset(particles[i]));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			CoreParticle* q = nodes[j]->key;
			if (Distance(p, q) <= thres)
			{
				merge(nodes[i], nodes[j]);
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

/*
Initialize a set of core clusters.
*/
vector<core_cluster*>
clusterCoreParticles(vector<CoreParticle*>& particles, float thres)
{
	vector<vector<CoreParticle*>> cores = clusterParticles(particles, thres);
	vector<core_cluster*> group(cores.size());
	CoreClusterFactory& factory = CoreClusterFactory::getInstance();
	for (int i = 0; i < cores.size(); ++i)
	{
		group[i] = factory.makeCoreCluster();
		group[i]->particles = cores[i];
	}
	return group;
}


/*
Check if there is a neighbor particle of the same generation that still has a room to grow. 
If so, then use the particle as a decendent of the core -> return true.
Else, this core is a legitimate one and should not be removed -> return false.
*/
bool
assignExit(vector<CoreParticle*> core, //pass a copy
			vector<CoreParticle*>& M, 
			int generation,
			int ndim,
			const int* dims)
{
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeEightNeighborhood(ndim);
	while (core.empty() == false)
	{
		vector<CoreParticle*> core2;
		for (int i = 0; i < core.size(); ++i)
		{
			vector<int> idx = coreParticle2Index(core[i], ndim);
			//int loc = Sub2Ind(idx, ndim, dims);
			float mind = std::numeric_limits<float>::infinity();
			CoreParticle* px = NULL;
			for (int n = 0; n < nbh.size(); ++n)
			{
				if (NeighborCheck(idx.begin(), nbh[n].begin(), ndim, dims))
				{
					CoreParticle* q = M[Sub2Ind(idx, nbh[n], ndim, dims)];
					if (q && q->gen == core[i]->gen && q->descendents.empty() == false)
					{
						float d = Length(idx, ndim);
						if (d < mind)
						{
							mind = d;
							px = q;
						}
					}
				}
			}
			if (px == NULL)
			{
				core2.push_back(core[i]);
			}
			else
			{
				if (px->ascendents.find(core[i]) == px->ascendents.end())
				{
					core[i]->descendents.insert(px);
					px->ascendents.insert(core[i]);
					px->trace_time[core[i]] = generation;
				}
			}
		}
		if (core.size() == core2.size())
		{
			break; //not possible to find exit for all particles.
		}
		core = core2;
	}
	return core.empty();
}

/*
Remove core particles (those with no descendents) that are too small (less than THRES).
*/
vector<CoreParticle*>
removeIsolatedCore(vector<CoreParticle*> core, //pass a copy
					int thres, vector<CoreParticle*>& M, 
					int generation,
					int ndim,
					const int* dims)
{
	vector<vector<CoreParticle*>> group = clusterParticles(core, 1.0f);
	//for each cluster, check if there is a neighbor that can be grown.
	//If so, assign its decentent to the neighbor (done all in assignExit).
	for (int g = 0; g < group.size(); ++g)
	{
		if (group[g].size() < thres)
		{
			if (assignExit(group[g], M, generation, ndim, dims) == true) //there is an exit particle
			{
				for (int k = 0; k < group[g].size(); ++k)
				{
					core.erase(find(core.begin(), core.end(), group[g][k]));
				}
			}
		}
	}
	return core;
}

/*
from each surface particle, move towrad the center and establish ascendent/descendent relation.
*/
vector<CoreParticle*>
propagateParticles(vector<int>& S,
					vector<unsigned char>& L,
					vector<CoreParticle*>& front,
					int mincore,
					int ndim,
					const int* dims)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> particles;
	vector<CoreParticle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	vector<CoreParticle*> Q;
	for (int i = 0; i < front.size(); ++i)
	{
		Q.push_back(front[i]);
		SetVoxel(mp, front[i], front[i], ndim, dims);
		particles.push_back(front[i]);
	}
	//int mincore = 5;
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeEightNeighborhood(ndim, dims);
	int gen = 1;
	while (Q.empty() == false)
	{
		vector<CoreParticle*> core;
		for (int i = 0; i < Q.size(); ++i)
		{
			SetVoxel(S, Q[i], gen, ndim, dims);
		}
		//G.push_back(Q); //separate per generation

		vector<CoreParticle*> Q2;
		for (int i = 0; i<Q.size(); ++i)
		{
			vector<int> sub = coreParticle2Index(Q[i], ndim);
			//int loc = Sub2Ind(coreParticle2Index(Q[i], ndim), ndim, dims);
			for (int n = 0; n < nbh.size(); ++n)
			{
				if (NeighborCheck(sub.begin(), nbh[n].begin(), ndim, dims))
				{
					int idx = Sub2Ind(sub, nbh[n], ndim, dims);
					if (L[idx] && S[idx] == 0)
					{
						CoreParticle* p = mp[idx];
						if (p == NULL)
						{
							p = coreParticleNdim(SubWithOffset(sub, nbh[n], ndim, dims), ndim);
							p->gen = gen;
							mp[idx] = p;
							Q2.push_back(p);
							particles.push_back(p);
						}
						if (p->ascendents.find(Q[i]) == p->ascendents.end())
						{
							p->ascendents.insert(Q[i]);
							Q[i]->descendents.insert(p);
							p->trace_time[Q[i]] = gen;
						}
					}
				}
			}
			if (Q[i]->descendents.empty())
			{
				core.push_back(Q[i]);
			}
		}
		core = removeIsolatedCore(core, mincore, mp, gen, ndim, dims);

		Q = Q2;
		gen++;
	}
	return particles;
}

/*
Erode an image (L) while not touching those in MASK.
*/
bool
erode(vector<unsigned char>& L, //foreground to be eroded.
		const vector<int>& mask, //non-negative voxels are untouched by erosion.
		int ndim,
		const int* dims)
{
	vector<CParticle4D> rem;
	if (ndim == 2)
	{
		for (int i = 0; i < dims[1]; ++i)
		{
			for (int j = 0; j < dims[0]; ++j)
			{
				unsigned char val = GetData2(L, j, i, dims[0], dims[1], (unsigned char)0);
				int val2 = GetData2(mask, j, i, dims[0], dims[1], -1); //check if it is an initial core
				if (val && val2 == -1)
				{
					if (GetData2(L, j - 1, i, dims[0], dims[1], (unsigned char)1) == 0 ||
						GetData2(L, j + 1, i, dims[0], dims[1], (unsigned char)1) == 0 ||
						GetData2(L, j, i - 1, dims[0], dims[1], (unsigned char)1) == 0 ||
						GetData2(L, j, i + 1, dims[0], dims[1], (unsigned char)1) == 0)
					{
						rem.push_back(CParticle4D(j, i));
					}

				}
			}
		}
	}
	else if (ndim == 3)
	{
		for (int k = 0; k < dims[2]; ++k)
		{
			for (int i = 0; i < dims[1]; ++i)
			{
				for (int j = 0; j < dims[0]; ++j)
				{
					unsigned char val = GetData3(L, j, i, k, dims[0], dims[1], dims[2], (unsigned char)0);
					int val2 = GetData3(mask, j, i, k, dims[0], dims[1], dims[2], -1); //check if it is an initial core
					if (val && val2 == -1)
					{
						if (
							GetData3(L, j - 1, i, k, dims[0], dims[1], dims[2], (unsigned char)1) == 0 ||
							GetData3(L, j + 1, i, k, dims[0], dims[1], dims[2], (unsigned char)1) == 0 ||
							GetData3(L, j, i - 1, k, dims[0], dims[1], dims[2], (unsigned char)1) == 0 ||
							GetData3(L, j, i + 1, k, dims[0], dims[1], dims[2], (unsigned char)1) == 0 ||
							GetData3(L, j, i, k - 1, dims[0], dims[1], dims[2], (unsigned char)1) == 0 ||
							GetData3(L, j, i, k + 1, dims[0], dims[1], dims[2], (unsigned char)1) == 0
							)
						{
							rem.push_back(CParticle4D(j, i, k));
						}

					}
				}
			}
		}
	}
	else if (ndim == 4)
	{
		for (int m = 0; m < dims[3]; ++m)
		{
			for (int k = 0; k < dims[2]; ++k)
			{
				for (int i = 0; i < dims[1]; ++i)
				{
					for (int j = 0; j < dims[0]; ++j)
					{
						unsigned char val = GetData4(L, j, i, k, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0);
						int val2 = GetData4(mask, j, i, k, m, dims[0], dims[1], dims[2], dims[3], -1); //check if it is an initial core
						if (val && val2 == -1)
						{
							if (
								GetData4(L, j - 1, i, k, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0 ||
								GetData4(L, j + 1, i, k, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0 ||
								GetData4(L, j, i - 1, k, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0 ||
								GetData4(L, j, i + 1, k, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0 ||
								GetData4(L, j, i, k - 1, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0 ||
								GetData4(L, j, i, k + 1, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0 ||
								GetData4(L, j, i, k, m - 1, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0 ||
								GetData4(L, j, i, k, m + 1, dims[0], dims[1], dims[2], dims[3], (unsigned char)1) == 0
								)
							{
								rem.push_back(CParticle4D(j, i, k, m));
							}

						}
					}
				}
			}
		}
	}
	else
	{
		mexErrMsgTxt("erode: unsupported number of dimensions. It has to be between 2 and 4.");
	}
	for (int i = 0; i < rem.size(); ++i)
	{
		SetVoxel(L, rem[i], (unsigned char)0, ndim, dims);
	}

	return !rem.empty();
}

/*
For each cluster, find the number of erosion operations to expose it.
*/
map<core_cluster*, int> 
coreWidth(vector<core_cluster*>& init, //initial groups of core particles
				const vector<int>& occupied, 
				int ndim,
				const int* dims)
{
	int nelm = numberOfElements(ndim, dims);
	vector<int> mask(nelm, -1);
	for (int i = 0; i < init.size(); ++i)
	{
		for (int j = 0; j < init[i]->size(); ++j)
		{
			SetVoxel(mask, init[i]->particles[j], i, ndim, dims); //used later in erosion
			init[i]->particles[j]->label = i;
		}
	}

	map<core_cluster*, int> width;
	vector<unsigned char> L(nelm, (unsigned char)0);
	for (int i = 0; i < occupied.size(); ++i)
	{
		if (occupied[i] >= 0)
		{
			L[i] = 1;
		}
	}
	int iter = 0;
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeEightNeighborhood(ndim, dims);
	while (true)
	{
		//perform erosion - but do not erode initial cores (stored in mask).
		bool bChanged = erode(L, mask, ndim, dims);
		if (!bChanged) break;
		iter++;

		for (int i = 0; i < init.size(); ++i)
		{
			if (width.find(init[i]) == width.end())
			{
				for (int j = 0; j < init[i]->size(); ++j) {
					CoreParticle* p = init[i]->particles[j];
					vector<int> sub = coreParticle2Index(p, ndim);
					if (onSurface(L, sub, nbh, ndim, dims))
					{
						width[init[i]] = iter;
						break;
					}
				}
			}
		}
	}

	return width;
}

/*
For each adjacent pair of clusters, find the cut cost, i.e. the number of erosion operations to separate the two.
*/
map<cluster_pair*, int>
clusterCutCost(vector<cluster_pair*>& pairs,
vector<core_cluster*>& init, //initial groups of core particles
vector<int>& occupied,
int ndim,
const int* dims)
{
	int nelm = numberOfElements(ndim, dims);
	vector<int> mask(nelm, -1);
	for (int i = 0; i < init.size(); ++i)
	{
		for (int j = 0; j < init[i]->size(); ++j)
		{
			SetVoxel(mask, init[i]->particles[j], i, ndim, dims); //used later in erosion
			init[i]->particles[j]->label = i;
		}
	}

	map<cluster_pair*, int> cuts;
	vector<unsigned char> L(nelm, (unsigned char)0);
	set<int> labels;
	for (int i = 0; i < occupied.size(); ++i)
	{
		if (occupied[i] >= 0)
		{
			L[i] = 1;
			labels.insert(occupied[i]);
		}
	}
	int iter = 0;
	vector<int> nbh = MakeFourNeighborhood(ndim, dims);
	while (true)
	{
		vector<int> C(L.size(), 0);
		int nc = ConnectedComponentAnalysisBigger(C, L, NeighborhoodFour, (unsigned char)0, ndim, dims);
		vector<Node<int>*> nodes;
		map<int, Node<int>*> nmap;
		for (set<int>::iterator it = labels.begin(); it != labels.end(); ++it)
		{
			Node<int>* n = makeset(*it);
			nodes.push_back(n);
			nmap[n->key] = n;
		}

		//for each cluster, check for the source (using occupied)
		vector<set<int>> group(nc + 1);
		for (int i = 0; i < C.size(); ++i)
		{
			int a = C[i];
			if (a > 0)
			{
				group[a].insert(occupied[i]);
			}
		}
		for (int i = 0; i < group.size(); ++i)
		{
			for (set<int>::iterator it = group[i].begin(); it != group[i].end(); ++it)
			{
				set<int>::iterator jt = it;
				jt++;
				for (; jt != group[i].end(); ++jt)
				{
					merge(nmap[*it], nmap[*jt]);
				}
			}
		}
		for (vector<cluster_pair*>::iterator it = pairs.begin(); it != pairs.end(); ++it)
		{
			cluster_pair* pr = *it;
			if (cuts.find(pr) == cuts.end())
			{
				if (findset(nmap[pr->u->id]) != findset(nmap[pr->v->id]))
				{
					cuts[pr] = iter;
				}
			}
		}
		for (int i = 0; i < nodes.size(); ++i)
		{
			delete nodes[i];
		}

		//perform erosion - but do not erode initial cores (stored in mask).
		bool bChanged = erode(L, mask, ndim, dims);
		if (!bChanged) break;
		iter++;
	}

	return cuts;
}

map<cluster_pair*, int>
clusterSeparationLength(const vector<CoreParticle*>& particles, //all particles
						const vector<core_cluster*>& init, //initial groups of core particles
						vector<int>& occupied, //mark the region growth and its origin
						int ndim,
						const int* dims)
{
	int nelm = numberOfElements(ndim, dims);
	vector<int> M(nelm, 0); //to store region-growth progress
	vector<CoreParticle*> pmap(nelm, NULL);
	for (int i = 0; i < particles.size(); ++i)
	{
		SetVoxel(pmap, particles[i], (CoreParticle*)particles[i], ndim, dims);
		particles[i]->label = -1; //label for non-core
	}
	vector<CoreParticle*> Q;
	map<int, core_cluster*> cmap;
	for (int i = 0; i < init.size(); ++i)
	{
		cmap[i] = init[i];
		for (int j = 0; j < init[i]->size(); ++j)
		{
			SetVoxel(occupied, init[i]->particles[j], i, ndim, dims);
			init[i]->particles[j]->label = init[i]->id;
			Q.push_back(init[i]->particles[j]);
		}
	}
	CoreClusterFactory& factory = CoreClusterFactory::getInstance();
	vector<cluster_pair*> pairs;
	map<cluster_pair*, int> separation; //maps adjacent pairs and its separation (based on the time it took for them to meet).
	vector<vector<int>> nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeFourNeighborhood(ndim, dims);
	int m = 0;
	while (Q.empty() == false)
	{
		m++;
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			SetVoxel(M, p, m, ndim, dims);
			vector<int> sub = coreParticle2Index(p, ndim);
			for (int n = 0; n < nbh.size(); ++n)
			{
				if (NeighborCheck(sub.begin(), nbh[n].begin(), ndim, dims))
				{
					int idx = Sub2Ind(sub, nbh[n], ndim, dims);
					CoreParticle* q = pmap[idx];
					if (q != NULL)
					{
						int k = GetVoxel(occupied, q, -1, ndim, dims);
						if (k < 0)
						{
							q->label = p->label;
							SetVoxel(occupied, q, p->label, ndim, dims);
							Q2.insert(q);
						}
						else if (k != p->label)
						{
							cluster_pair* pr = factory.makeClusterPair(cmap[p->label], cmap[q->label]);
							if (find(pairs.begin(), pairs.end(), pr) == pairs.end())
							{
								pairs.push_back(pr);
								separation[pr] = m;
							}
						}
					}
				}
			}
		}
		Q.clear();
		Q.insert(Q.begin(), Q2.begin(), Q2.end());
	}
	return separation;
}
/*
This routine first establish a graph by region-growing from each core - when two regions originated from different cores touch,
then the cores are considered adjacent.
It then perform grouping using the graph -- if the edge strength is less than thres, then the nodes are put into the same group.
*/
vector<core_cluster*>
groupCoreParticles(vector<CoreParticle*>& particles, //all particles
					vector<core_cluster*>& init, //initial groups of core particles
					//int nclusters, //desired number of clusters 
					float thres, //cut off point
					int  ndim,
					const int* dims)
{
	vector<int> occupied(numberOfElements(ndim, dims), -1);
	map<cluster_pair*, int> separation = clusterSeparationLength(particles, init, occupied, ndim, dims);
	vector<cluster_pair*> core_pairs; //pairs of adjacent core clusters
	for (map<cluster_pair*, int>::iterator it = separation.begin(); it != separation.end(); ++it)
	{
		core_pairs.push_back(it->first);
	}
	map<cluster_pair*, int> cuts = clusterCutCost(core_pairs, init, occupied, ndim, dims);

	//This should be done when clusters are first formed...
	map<core_cluster*, int> width = coreWidth(init, occupied, ndim, dims);
	for (int i = 0; i<init.size(); ++i)
	{
		init[i]->width = width[init[i]];
	}

	vector<core_cluster*> clusters = init;
	int npairs = core_pairs.size();
	CoreClusterFactory& factory = CoreClusterFactory::getInstance();
	for (int i = 0; i < npairs; ++i)
	{
		cluster_pair* p = core_pairs[i];
		p->length = separation[p];
		p->cut = cuts[p];
		p->u->pairs.push_back(p);
		//make a copy with opposite direction
		cluster_pair* q = factory.makeClusterPair(p->v, p->u);
		q->length = p->length;
		q->cut = p->cut;
		q->u->pairs.push_back(q);

		CParticle4D cu = q->u->centroid();
		CParticle4D cv = q->v->centroid();
		/*printf("Pair: %d (%d:%d,%d) and %d (%d:%d,%d) - %d, %d, %d, %d, %f\n",
			q->u->id, q->u->size(), cu.m_X + 1, cu.m_Y + 1, q->v->id, q->v->size(), cv.m_X + 1, cv.m_Y + 1, (int)q->length, (int)q->cut,
			q->u->width, q->v->width, q->weight());*/
	}

	while (clusters.size() > 0) //nclusters)
	{
		pair<int, int> choice;
		float mincost = std::numeric_limits<float>::infinity();
		for (int i = 0; i < clusters.size(); ++i)
		{
			for (int j = 0; j < clusters[i]->pairs.size(); ++j)
			{
				if (clusters[i]->pairs[j]->weight() < mincost)
				{
					mincost = clusters[i]->pairs[j]->weight();
					choice = pair<int, int>(i, j);
				}
			}
		}
		if (mincost <= thres) //std::numeric_limits<float>::infinity())
		{
			cluster_pair* pr = clusters[choice.first]->pairs[choice.second];
			CParticle4D cu = pr->u->centroid();
			CParticle4D cv = pr->v->centroid();
			printf("Merging: %d (%d,%d) and %d (%d,%d) - %d, %d, %d, %d, %f\n",
				pr->u->id, cu.m_X + 1, cu.m_Y + 1, pr->v->id, cv.m_X + 1, cv.m_Y + 1, (int)pr->length, (int)pr->cut,
				pr->u->width, pr->v->width, pr->weight());
			core_cluster* c = merge_cores(pr->u, pr->v);
			clusters.erase(find(clusters.begin(), clusters.end(), pr->u));
			clusters.erase(find(clusters.begin(), clusters.end(), pr->v));
			clusters.push_back(c);
		}
		else
		{
			break;
		}
	}

	//reset label
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
	}
	return clusters;
}

set<CoreParticle*> trace(CoreParticle* p,
						map<CoreParticle*,int>& arrival)
{
	set<CoreParticle*> S;
	set<CoreParticle*> Q;
	Q.insert(p);
	while (Q.empty() == false)
	{
		set<CoreParticle*> next;
		for (set<CoreParticle*>::iterator it = Q.begin(); it != Q.end(); ++it)
		{
			CoreParticle* q = *it;
			if (q->label == 0) q->label = -2;
			S.insert(q);
			for (set<CoreParticle*>::iterator it = q->ascendents.begin(); it != q->ascendents.end(); ++it)
			{
				CoreParticle* a = *it;
				next.insert(a);
				if (arrival.find(a) == arrival.end())
				{
					a->label = q->label;
					arrival[a] = q->trace_time[a];
				}
				else if (arrival[a] > q->trace_time[a])
				{
					a->label = q->label;
					arrival[a] = q->trace_time[a];
				}
				else if (arrival[a] == q->trace_time[a] && a->label != q->label)
				{
					a->label = 0;
				}
				//else if(a->label != q->label)  a->label = -1; //conflict.
			}
		}
		Q = next;
	}
	return S;
}

set<CoreParticle*>
colorParticles(vector<core_cluster*>& cores)
{
	set<CoreParticle*> S;
	map<CoreParticle*, int> arrival;
	int idx = 0;
	int count = 0;
	for (int i = 0; i < cores.size(); ++i)
	{
		for (int j = 0; j < cores[i]->size(); ++j)
		{
			CoreParticle* p = cores[i]->particles[j];
			arrival[p] = 0;
		}
		if (count < cores[i]->size())
		{
			idx = i;
			count = cores[i]->size();
		}
	}
	for (int i = 0; i < cores.size(); ++i)
	{
		//int i = idx;
		int label = i + 1;
		//printf("There are %d particles in core #%d.\n", cores[i]->size(), i);
		for (int j = 0; j < cores[i]->size(); ++j)
		{
			CoreParticle* p = cores[i]->particles[j];
			p->label = label;
			set<CoreParticle*> s = trace(p, arrival);
			S.insert(s.begin(), s.end());
			//printf("%d(%d,%d,%d)\n", p->id, p->x, p->y, p->z);
		}
	}
	return S;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{   
	if (nrhs < 1|| nlhs < 0)
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
	float costCutoff = 10.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(costCutoff, prhs[1], classMode);
	}
	int minsize = 5;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minsize, prhs[2], classMode);
	}

	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	int nvoxels = numberOfElements(ndimL, dimsL);

	vector<int> S(nvoxels, 0);

	vector<CoreParticle*> front = ExtactSurface(L, ndimL, dimsL);
	vector<CoreParticle*> particles = propagateParticles(S, L, front, minsize, ndimL, dimsL);
	vector<CoreParticle*> core = pickCore(particles);
	vector<core_cluster*> init = clusterCoreParticles(core, 1.1); //initial clusters
	vector<core_cluster*> coregroup = groupCoreParticles(particles, init, costCutoff, ndimL, dimsL);
	set<CoreParticle*> colored = colorParticles(coregroup);

	if (nlhs >= 1)
	{
		int count = 0;
		for (int i = 0; i < init.size(); ++i)
		{
			count += init[i]->size();
		}
		int dims[] = { count, 9 };
		map<CoreParticle*, int> cmap;
		for (int i = 0, k = 0; i < coregroup.size(); ++i)
		{
			for (int j = 0; j < coregroup[i]->size(); ++j, ++k)
			{
				CoreParticle* p = coregroup[i]->particles[j];
				cmap[p] = i;
			}
		}

		vector<int> F(dims[0] * dims[1]);
		for (int i = 0, k = 0; i < init.size(); ++i)
		{
			for (int j = 0; j < init[i]->size(); ++j, ++k)
			{
				CoreParticle* p = init[i]->particles[j];
				SetData2(F, k, 0, dims[0], dims[1], p->x);
				SetData2(F, k, 1, dims[0], dims[1], p->y);
				SetData2(F, k, 2, dims[0], dims[1], p->z);
				SetData2(F, k, 3, dims[0], dims[1], p->t);
				SetData2(F, k, 4, dims[0], dims[1], p->label);
				SetData2(F, k, 5, dims[0], dims[1], p->gen);
				SetData2(F, k, 6, dims[0], dims[1], p->id);
				SetData2(F, k, 7, dims[0], dims[1], i);
				SetData2(F, k, 8, dims[0], dims[1], cmap[p]);
			}
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { colored.size(), 7 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (set<CoreParticle*>::iterator it = colored.begin(); it != colored.end(); ++it, ++i)
		{
			CoreParticle* p = *it;
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
		plhs[2] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 4)
	{
		vector<int> T(numberOfElements(ndimL, dimsL), 0);
		for (set<CoreParticle*>::iterator it = colored.begin(); it != colored.end(); ++it)
		{
			CoreParticle* p = *it;
			SetVoxel(T, p, p->label + 1, ndimL, dimsL);
		}
		plhs[3] = StoreData(T, mxINT32_CLASS, ndimL, dimsL);

	}
	mexUnlock();
	factory.clean();
	CoreClusterFactory::getInstance().clean();
}

