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

struct particle_cluster;

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
		cluster = NULL;
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
	particle_cluster* cluster;
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

struct cluster_pair;
struct particle_cluster
{
	particle_cluster(vector<CoreParticle*>& vp, int gen0, int id0)
	{
		particles = vp;
		gen = gen0;
		id = id0;
		width = 0;
		label = 0;
		value = 0.f;
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
	set<particle_cluster*> ascendents;
	set<particle_cluster*> descendents;
	int size() { return particles.size(); }
	int width;
	int gen; //generation
	int id;
	int label; //for general purpose
	float value; //for general purpose
};

struct cluster_pair
{
	cluster_pair() {
		u = NULL;
		v = NULL;
		cut = std::numeric_limits<float>::infinity();
		length = 0;
	}
	particle_cluster* u;
	particle_cluster* v;
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
	particle_cluster* makeCoreCluster(vector<CoreParticle*>& vp, int gen)
	{
		particle_cluster* c = new particle_cluster(vp, gen, _id++);
		clusters.push_back(c);
		return c;
	}
	cluster_pair* makeClusterPair(particle_cluster* a, particle_cluster* b)
	{
		cluster_pair* p = new cluster_pair();
		p->u = a;
		p->v = b;
		pairs.push_back(p);
		return p;
	}
	void clean()
	{
		for (int i = 0; i < clusters.size(); ++i)
		{
			delete clusters[i];
		}
		clusters.clear();
		for (int i = 0; i < pairs.size(); ++i)
		{
			delete pairs[i];
		}
		pairs.clear();
		_id = 0;
	}
	vector<particle_cluster*> clusters;
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

/*particle_cluster*
merge_cores(particle_cluster* a, particle_cluster* b)
{
	CoreClusterFactory& factory = CoreClusterFactory::getInstance();
	particle_cluster* c = factory.makeCoreCluster();
	c->particles = a->particles;
	c->particles.insert(c->particles.end(), b->particles.begin(), b->particles.end());
	c->width = Max(a->width, b->width);  //Is this right??
	for (int i = 0; i < a->pairs.size(); ++i)
	{
		if (a->pairs[i]->v != b)
		{
			a->pairs[i]->u = c;
			particle_cluster* v = a->pairs[i]->v;
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
			particle_cluster* v = b->pairs[i]->v;
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
}*/

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
Length(float x, float y, float z = 0, float t = 0)
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

set<particle_cluster*>
collectAscendents(particle_cluster* c)
{
	//collect all ascendents of pc including itself.
	set<particle_cluster*> all;
	set<particle_cluster*> next;
	next.insert(c);
	while (next.empty() == false)
	{
		set<particle_cluster*> Q;
		for (set<particle_cluster*>::iterator it = next.begin(); it != next.end(); ++it)
		{
			all.insert(*it);
			for (set<particle_cluster*>::iterator jt = (*it)->ascendents.begin(); jt != (*it)->ascendents.end(); ++jt)
			{
				Q.insert(*it);

			}
		}
		next = Q;
	}
	return all;
}


bool
noCommonAscendents(particle_cluster* a, particle_cluster* b)
{
	set<particle_cluster*> sa = collectAscendents(a);
	set<particle_cluster*> sb = collectAscendents(b);
	vector<particle_cluster*> diff;
	set_difference(sa.begin(), sa.end(), sb.begin(), sb.end(), std::inserter(diff, diff.begin()));
	return diff.empty();
}

int 
numDisjointClusters(set<particle_cluster*>& P)
{
	set<CoreParticle*> S;
	for (set<particle_cluster*>::iterator it = P.begin(); it != P.end(); ++it)
	{
		S.insert((*it)->particles.begin(), (*it)->particles.end());
	}
	vector<Node<CoreParticle*>*> nodes;
	for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
	{
		nodes.push_back(makeset(*it));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		for (int j =i+1; j < nodes.size(); ++j)
		{
			CoreParticle* q = nodes[j]->key;
			if (Abs(p->x - q->x) <= 1 && Abs(p->y - q->y) <= 1 && Abs(p->z - q->z) <= 1 && Abs(p->t - q->t) <= 1)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<CoreParticle*>*> reps = clusters(nodes);
	int nc = reps.size();
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return nc;
}

bool
isAdjacent(CoreParticle* p, vector<CoreParticle*>& particles)
{
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* q = particles[i];
		int cnt = 0;
		cnt += Abs(p->x - q->x);
		cnt += Abs(p->y - q->y);
		cnt += Abs(p->z - q->z);
		cnt += Abs(p->t - q->t);
		if (cnt == 1)
		{
			return true;
		}
	}
	return false;
}

vector<particle_cluster*>
pickCore(vector<particle_cluster*>& pclusters, float thres)
{
	vector<particle_cluster*> core;
	for (int i = 0; i < pclusters.size(); ++i)
	{
		//set of ascendents
		set<CoreParticle*> aset;
		for (set<particle_cluster*>::iterator it = pclusters[i]->ascendents.begin(); it != pclusters[i]->ascendents.end(); ++it)
		{
			particle_cluster* pc = *it;
			for (int j = 0; j < pc->particles.size(); ++j) {
				if (isAdjacent(pc->particles[j], pclusters[i]->particles))
				{
					aset.insert(pc->particles[j]);
				}
			}
		}
		//set of descendents
		set<CoreParticle*> dset;
		for (set<particle_cluster*>::iterator it = pclusters[i]->descendents.begin(); it != pclusters[i]->descendents.end(); ++it)
		{
			particle_cluster* pc = *it;
			for (int j = 0; j < pc->particles.size(); ++j) {
				if (isAdjacent(pc->particles[j], pclusters[i]->particles))
				{
					dset.insert(pc->particles[j]);
				}
			}
		}
		float measure = (float)dset.size() / (float)aset.size();
		pclusters[i]->value = measure;
		//printf("corenesss: %d %d %d %f\n", i, dset.size(), aset.size(), measure);
		if (aset.size() > 0 && measure < thres)
		{
			core.push_back(pclusters[i]);
			pclusters[i]->label = 1;
		}
		else {
			pclusters[i]->label = 0;
		}
	}
	vector<particle_cluster*> core2;
	for (int i = 0; i < core.size(); ++i)
	{
		bool bCore = true;
		for (set<particle_cluster*>::iterator it = core[i]->descendents.begin(); it != core[i]->descendents.end(); ++it)
		{
			if ((*it)->label == 1)
			{
				bCore = false;
				break;
			}
		}
		if (bCore)
		{
			core2.push_back(core[i]);
		}
	}

	printf("Cores reduced from %d to %d clusters.\n", core.size(), core2.size());
	return core2;
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
from each surface particle, move towrad the center and establish ascendent/descendent relation.
*/
vector<particle_cluster*>
propagateParticles(vector<int>& S,
vector<unsigned char>& L,
vector<CoreParticle*>& front,
int mincore,
int ndim,
const int* dims)
{
	CoreClusterFactory& cfactory = CoreClusterFactory::getInstance();
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> particles;
	vector<particle_cluster*> pclusters;
	vector<CoreParticle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	vector<CoreParticle*> Q;
	for (int i = 0; i < front.size(); ++i)
	{
		Q.push_back(front[i]);
		SetVoxel(mp, front[i], front[i], ndim, dims);
		particles.push_back(front[i]);
	}
	vector<unsigned char> F(numberOfElements(ndim, dims)); //for connected components
	vector<int> C(numberOfElements(ndim, dims)); //for connected component labels
	//int mincore = 5;
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeEightNeighborhood(ndim, dims);
	int gen = 1;
	while (Q.empty() == false)
	{
		for (int i = 0; i < F.size(); ++i)
		{
			F[i] = 0;
			C[i] = 0;
		}
		vector<CoreParticle*> core;
		for (int i = 0; i < Q.size(); ++i)
		{
			SetVoxel(S, Q[i], gen, ndim, dims);
			SetVoxel(F, Q[i], (unsigned char)1, ndim, dims);
		}
		int nc = ConnectedComponentAnalysisBigger(C, F, NeighborhoodFour, (unsigned char)0, ndim, dims);
		vector<vector<CoreParticle*>> vvp(nc);
		for (int i = 0; i < C.size(); ++i)
		{
			if (C[i] > 0)
			{
				CoreParticle* p = mp[i];
				vvp[C[i] - 1].push_back(p);
			}
		}
		vector<particle_cluster*> vpc;
		for (int i = 0; i < nc; ++i)
		{
			if (vvp[i].empty()) continue;

			particle_cluster* pc = cfactory.makeCoreCluster(vvp[i],  gen);
			for (int j = 0; j < vvp[i].size(); ++j)
			{
				vvp[i][j]->cluster = pc;
			}
			vpc.push_back(pc);
		}
		for (int i = 0; i < vpc.size(); ++i)
		{
			set<particle_cluster*> spc;
			for (int j = 0; j < vpc[i]->particles.size(); ++j)
			{
				CoreParticle* p = vpc[i]->particles[j];
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					spc.insert(q->cluster);
				}
			}
			for (set<particle_cluster*>::iterator jt = spc.begin(); jt != spc.end(); ++jt)
			{
				vpc[i]->ascendents.insert(*jt);
				(*jt)->descendents.insert(vpc[i]);
			}
		}
		pclusters.insert(pclusters.end(), vpc.begin(), vpc.end());

		vector<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			vector<int> sub = coreParticle2Index(Q[i], ndim);
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
						}
					}
				}
			}
			/*if (Q[i]->descendents.empty())
			{
			core.push_back(Q[i]);
			}*/
		}
		//core = removeIsolatedCore(core, mincore, mp, gen, ndim, dims);

		Q = Q2;
		gen++;
	}
	return pclusters;
}

void
colorParticles(vector<particle_cluster*>& clusters, vector<particle_cluster*>& cores)
{
	for (int i = 0; i < clusters.size(); ++i)
	{
		clusters[i]->label = 0;
	}
	set<particle_cluster*> Q;
	for (int i = 0; i < cores.size(); ++i)
	{
		cores[i]->label = i+1;
		Q.insert(cores[i]);
		printf("Core %d: %d, %d\n", i + 1, cores[i]->gen, cores[i]->size());
	}
	//Q.insert(cores[cores.size() - 1]);

	while (Q.empty()==false)
	{
		set<particle_cluster*> Q2;
		for (set<particle_cluster*>::iterator it = Q.begin(); it != Q.end(); ++it)
		{
			particle_cluster* pc = *it;
			for (set<particle_cluster*>::iterator jt = pc->ascendents.begin(); jt != pc->ascendents.end(); ++jt)
			{
				particle_cluster* pd = *jt;
				if (pd->label == 0)
				{
					pd->label = pc->label;
					Q2.insert(pd);
				}
				/*else if (pc->gen > pd->gen)
				{
					pd->label = pc->label;
					Q2.insert(pd);
				}*/
			}
		}
		Q = Q2;
	}
	for (int i = 0; i < clusters.size(); ++i)
	{
		for (int j = 0; j < clusters[i]->particles.size(); ++j)
		{
			clusters[i]->particles[j]->label = clusters[i]->label;
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
	float costCutoff = 0.3f;
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
	vector<particle_cluster*> clusters = propagateParticles(S, L, front, minsize, ndimL, dimsL);
	vector<particle_cluster*> core = pickCore(clusters, costCutoff);
	colorParticles(clusters, core);
	vector<int> T(nvoxels, 0);
	for (int i = 0; i < clusters.size(); ++i)
	{
		for (int j = 0; j < clusters[i]->particles.size(); ++j)
		{
			CoreParticle* p = clusters[i]->particles[j];
			SetVoxel(T, p, p->label, ndimL, dimsL);
		}
	}

	if (nlhs >= 1)
	{
		int count = 0;
		for (int i = 0; i < core.size(); ++i)
		{
			count += core[i]->size();
		}
		int dims[] = { count, 8 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0, k = 0; i < core.size(); ++i)
		{
			for (int j = 0; j < core[i]->size(); ++j, ++k)
			{
				CoreParticle* p = core[i]->particles[j];
				SetData2(F, k, 0, dims[0], dims[1], p->x);
				SetData2(F, k, 1, dims[0], dims[1], p->y);
				SetData2(F, k, 2, dims[0], dims[1], p->z);
				SetData2(F, k, 3, dims[0], dims[1], p->t);
				SetData2(F, k, 4, dims[0], dims[1], p->label);
				SetData2(F, k, 5, dims[0], dims[1], p->gen);
				SetData2(F, k, 6, dims[0], dims[1], p->id);
				SetData2(F, k, 7, dims[0], dims[1], p->cluster->id);
			}
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int count = 0;
		for (int i = 0; i < clusters.size(); ++i)
		{
			count += clusters[i]->size();
		}
		int dims[] = { count, 8 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0, k = 0; i < clusters.size(); ++i)
		{
			for (int j = 0; j < clusters[i]->size(); ++j, ++k)
			{
				CoreParticle* p = clusters[i]->particles[j];
				SetData2(F, k, 0, dims[0], dims[1], p->x);
				SetData2(F, k, 1, dims[0], dims[1], p->y);
				SetData2(F, k, 2, dims[0], dims[1], p->z);
				SetData2(F, k, 3, dims[0], dims[1], p->t);
				SetData2(F, k, 4, dims[0], dims[1], p->label);
				SetData2(F, k, 5, dims[0], dims[1], p->gen);
				SetData2(F, k, 6, dims[0], dims[1], p->id);
				SetData2(F, k, 7, dims[0], dims[1], p->cluster->id);
			}
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 3)
	{
		plhs[2] = StoreData(S, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 4)
	{
		plhs[3] = StoreData(T, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 5)
	{
		int dims[] = {core.size(), 1 };
		vector<vector<int>> F;
		for (int i = 0, k = 0; i < core.size(); ++i)
		{
			vector<int> vids;
			set<particle_cluster*> Q;
			Q.insert(core[i]);
			while (Q.empty() == false)
			{
				set<particle_cluster*> Q2;
				for (set<particle_cluster*>::iterator it = Q.begin(); it != Q.end(); ++it)
				{
					particle_cluster* pc = *it;
					vids.push_back(pc->id);
					for (set<particle_cluster*>::iterator jt = pc->ascendents.begin(); jt != pc->ascendents.end(); ++jt)
					{
						Q2.insert(*jt);
					}
				}
				Q = Q2;
			}
			F.push_back(vids);
		}
		plhs[4] = StoreDataCell(F, mxINT32_CLASS, 2, dims, 1);
	}
	mexUnlock();
	factory.clean();
	CoreClusterFactory::getInstance().clean();
}

