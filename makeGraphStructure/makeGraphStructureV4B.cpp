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
#include <Kruskal.h>
#include <Dijkstra.h>

bool _EightNeighbor = false;
bool _Use_New = true;
int _debug_id = 302; 
int _debug_id2 = 303;
float strong_core_thres = 0.25;

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
		dval = 0;
		src = NULL;
		//selected = false;
		core = -1; //this will be set either 0, 1, or 2.
		vertex = NULL;
	}
	int x;
	int y;
	int z;
	int t;
	int label;
	int gen; //generation
	int id;
	float value; //generic value
	float dval; //distance value
	//bool selected;
	int core;

	set<CoreParticle*> ascendents;
	set<CoreParticle*> descendents;
	vector<CoreParticle*> neighbors;
	vector<CoreParticle*> neighbors8;
	//CoreParticle* pi; //path from the core
	CoreParticle* src; //strong core of this particle
	Vertex<CoreParticle*>* vertex;
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

float
length(float x, float y, float z, float t)
{
	return sqrt(x*x + y*y + z*z + t*t);
}

vector<vector<int>> MakeEightNeighborhoodHere(int n)
{
	if (n == 2)
	{
		vector<vector<int>> idx(8);
		int xoff[] = { 0, -1, 1, 0, -1, 1, -1, 1 };
		int yoff[] = { -1, 0, 0, 1, -1, -1, 1, 1 };
		for (int i = 0; i < idx.size(); ++i)
		{
			vector<int> jdx(2);
			jdx[0] = xoff[i];
			jdx[1] = yoff[i];
			idx[i] = jdx;
		}
		return idx;
	}
	else if (n == 3)
	{
		vector<pair<float,CParticle4D>> vals;
		for (int x = -1; x <= 1; ++x)
		{
			for (int y = -1; y <= 1; ++y)
			{
				for (int z = -1; z <= 1; ++z)
				{
					if (!(x == 0 && y == 0 && z == 0))
					{
						CParticle4D p(x, y, z, 0);
						float len = length(x, y, z, 0);
						vals.push_back(pair<float, CParticle4D>(len, p));
					}
				}
			}
		}
		sort(vals.begin(), vals.end());
		vector<vector<int>> idx;
		for (int i = 0; i < vals.size(); ++i)
		{
			vector<int> jdx(3);
			jdx[0] = vals[i].second.m_X;
			jdx[1] = vals[i].second.m_Y;
			jdx[2] = vals[i].second.m_Z;
			idx.push_back(jdx);
		}
		return idx;
	}
	else
	{
		return MakeEightNeighborhood(n);
	}
}

struct NeighborhoodFactory
{
public:
	static NeighborhoodFactory& getInstance(int n = 0)
	{
		static NeighborhoodFactory instance(n);
		if (n > 0 && n != instance.ndim)
		{
			instance.neighbor4 = MakeFourNeighborhood(n);
			instance.neighbor8 = MakeEightNeighborhoodHere(n);
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
	}
	~NeighborhoodFactory()
	{
	}
	NeighborhoodFactory(NeighborhoodFactory& f){}
	NeighborhoodFactory operator=(NeighborhoodFactory& f){}
};

struct BoundingBox
{
	BoundingBox(CParticle4D a, CParticle4D b)
	{
		minc = CParticle4D(Min(a.m_X, b.m_X), Min(a.m_Y, b.m_Y), Min(a.m_Z, b.m_Z), Min(a.m_T, b.m_T));
		maxc = CParticle4D(Max(a.m_X, b.m_X), Max(a.m_Y, b.m_Y), Max(a.m_Z, b.m_Z), Max(a.m_T, b.m_T));
	}
	BoundingBox(set<CoreParticle*>& sp)
	{
		int maxval = std::numeric_limits<int>::max();
		int minval = std::numeric_limits<int>::min();
		minc = CParticle4D(maxval, maxval, maxval, maxval);
		maxc = CParticle4D(minval, minval, minval, minval);
		for (set<CoreParticle*>::iterator it = sp.begin(); it != sp.end(); ++it)
		{
			minc.m_X = Min(minc.m_X, (*it)->x);
			minc.m_Y = Min(minc.m_Y, (*it)->y);
			minc.m_Z = Min(minc.m_Z, (*it)->z);
			minc.m_T = Min(minc.m_T, (*it)->t);
			maxc.m_X = Max(maxc.m_X, (*it)->x);
			maxc.m_Y = Max(maxc.m_Y, (*it)->y);
			maxc.m_Z = Max(maxc.m_Z, (*it)->z);
			maxc.m_T = Max(maxc.m_T, (*it)->t);
		}
	}
	float volume()
	{
		return (float)(maxc.m_X - minc.m_X + 1)*(maxc.m_Y - minc.m_Y + 1)*(maxc.m_Z - minc.m_Z + 1)*(maxc.m_T - minc.m_T + 1);
	}
	bool inside(CParticle4D& p)
	{
		return p.m_X >= minc.m_X && p.m_X <= maxc.m_X && p.m_Y >= minc.m_Y && p.m_Y <= maxc.m_Y &&
			p.m_Z >= minc.m_Z && p.m_Z <= maxc.m_Z && p.m_T >= minc.m_T && p.m_T <= maxc.m_T;
	}
	CParticle4D minc;
	CParticle4D maxc;
};

float overlapVolume(BoundingBox& b1, BoundingBox& b2)
{
	int x1 = Max(b1.minc.m_X, b2.minc.m_X);
	int x2 = Min(b1.maxc.m_X, b2.maxc.m_X);
	int y1 = Max(b1.minc.m_Y, b2.minc.m_Y);
	int y2 = Min(b1.maxc.m_Y, b2.maxc.m_Y);
	int z1 = Max(b1.minc.m_Z, b2.minc.m_Z);
	int z2 = Min(b1.maxc.m_Z, b2.maxc.m_Z);
	int t1 = Max(b1.minc.m_T, b2.minc.m_T);
	int t2 = Min(b1.maxc.m_T, b2.maxc.m_T);
	if (x1 > x2 || y1 > y2 || z1 > z2 || t1 > t2) return 0.0f;
	else {
		return (x2 - x1 + 1)*(y2 - y1 + 1)*(z2 - z1 + 1)*(t2 - t1 + 1);
	}
	/*CParticle4D minc, maxc;
	if (b1.inside(b2.minc)) minc = b2.minc;
	else  if (b2.inside(b1.minc)) minc = b1.minc;
	else return 0.0f; //no overlap
	if (b1.inside(b2.maxc)) maxc = b2.maxc;
	else if (b2.inside(b1.maxc)) maxc = b1.maxc;
	else return 0.0f; //no overlap

	BoundingBox b3(minc, maxc);
	return b3.volume();*/
}

float
dotProduct(float x1, float y1, float z1, float t1, float x2, float y2, float z2, float t2)
{
	return x1*x2 + y1*y2 + z1*z2 + t1*t2;
}

/*
When two particles (p and q) meet originated from different strong core sduring the growth, and a graph edge will be formed.
This routine provides the weight of the edge.
*/
float weight(CoreParticle* p, CoreParticle* q)
{
	//return sqrt(p->src->dval * q->src->dval) / sqrt(p->dval * q->dval);
	//float len = length(p->src->x - q->src->x, p->src->y - q->src->y, p->src->z - q->src->z, p->src->t - q->src->t);
	//return (len + sqrt(p->src->dval * q->src->dval) + 1.0) / sqrt(p->dval * q->dval);
	return length(p->src->x - q->src->x, p->src->y - q->src->y, p->src->z - q->src->z, p->src->t - q->src->t);
}

bool
surfaceParticle(CoreParticle* p, int ndim)
{
	int ns = NeighborhoodFactory::getInstance(ndim).neighbor4.size();
	if (_EightNeighbor)
	{
		ns = NeighborhoodFactory::getInstance(ndim).neighbor8.size();
	}
	return p->neighbors.size() < ns;
}

bool
medialParticle(CoreParticle* p, int ndim)
{
	return p->ascendents.size() >= 2 * (ndim - 1) || p->descendents.size() == 0;
}

vector<vector<CoreParticle*>>
clusterParticles(vector<CoreParticle*>& particles);

/*
Ascendents of p are not linearly separable.
*/
bool
strongMedialParticle(CoreParticle* p, vector<vector<CoreParticle*>>& groups, int ndim)
{
	if (p->core >= 0) return p->core > 0;
	if (p->descendents.size() == 0)
	{
		p->core = 2;
		//return true;
	}
	else
	{
		p->core = 0;
	}
	if (ndim == 2 || ndim == 3)
	{
		/*if (p->ascendents.size() >= 2 * (ndim - 1))
		{
			p->core = 1;
		}*/
		/*if (p->dval >= 7.0f)
		{
			bool blmax = true;
			for (int i = 0; i<p->neighbors8.size(); ++i)
			{
				if (p->neighbors8[i]->dval >= p->dval)
				{
					blmax = false;
					break;
				}
			}
			if (blmax) p->core = 1;
		}*/
	}
	return p->core > 0;
}

void
labelStrongCoreParticles(vector<CoreParticle*>& particles, int ndim)
{
	vector<CoreParticle*> surfaces;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (surfaceParticle(particles[i], ndim))
		{
			surfaces.push_back(particles[i]);
		}
	}
	vector<vector<CoreParticle*>> groups = clusterParticles(surfaces);
	map<CoreParticle*, int> imap;
	for (int i = 0; i < particles.size(); ++i)
	{
		strongMedialParticle(particles[i], groups, ndim);
	}
}

bool
inflectionParticle(CoreParticle* p, int ndim)
{
	return surfaceParticle(p, ndim) && p->descendents.size() >= 2;  //ndim;
}//


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
	NeighborhoodFactory& nfactory = NeighborhoodFactory::getInstance(ndim);
	vector<vector<int>> nbh = nfactory.neighbor4;
	vector<vector<int>> nbh8 = nfactory.neighbor8;

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
		for (int n = 0; n < nbh8.size(); ++n)
		{
			if (NeighborCheck(sub.begin(), nbh8[n].begin(), ndim, dims))
			{
				int idx = Sub2Ind(sub, nbh8[n], ndim, dims);
				if (mp[idx] != NULL)
				{
					p->neighbors8.push_back(mp[idx]);
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
		for (int j = 0; j < p->neighbors8.size(); ++j)
		{
			CoreParticle* q = p->neighbors8[j];
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


/*
Based on the L2 distance value, find the neighbor that is most straight.

This version uses 8-neighbors.
*/
CoreParticle*
representativeNeighbor(CoreParticle* p, int gen)
{
	float minerr = std::numeric_limits<float>::infinity();
	CoreParticle* best = NULL;
	for (int i = 0; i < p->neighbors8.size(); ++i)
	{
		CoreParticle* q = p->neighbors8[i];
		float len = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
		if (q->gen == gen)
		{
			float err = Abs(p->dval + len - q->dval);
			if (err < minerr)
			{
				minerr = err;
				best = q;
			}
		}
	}
	return best;
}

vector<CoreParticle*>
selectNeighbors(CoreParticle*p, int gen, bool bEight = true)
{
	vector<CoreParticle*> vp;
	if (bEight)
	{
		for (int i = 0; i < p->neighbors8.size(); ++i)
		{
			if (p->neighbors8[i]->gen == gen)
			{
				vp.push_back(p->neighbors8[i]);
			}
		}
	}
	else
	{
		for (int i = 0; i < p->neighbors.size(); ++i)
		{
			if (p->neighbors[i]->gen == gen)
			{
				vp.push_back(p->neighbors[i]);
			}
		}
	}
	return vp;
}

vector<float>
centroid(vector<CoreParticle*>& vp)
{
	vector<float> v(4, 0.0f);
	for (int j = 0; j < vp.size(); ++j)
	{
		v[0] += vp[j]->x;
		v[1] += vp[j]->y;
		v[2] += vp[j]->z;
		v[3] += vp[j]->t;
	}
	v[0] /= vp.size(); v[1] /= vp.size(); v[2] /= vp.size(); v[3] /= vp.size();
	return v;
}

/*
find a particle among VP that is closest to C with 4 elements.
*/
CoreParticle*
closestParticle(vector<CoreParticle*>& vp, vector<float> c)
{
	//vector<float> v = centroid(vp);
	CoreParticle* rep = NULL;
	float mind = std::numeric_limits<float>::infinity();
	for (int j = 0; j < vp.size(); ++j)
	{
		float d = length(c[0] - vp[j]->x, c[1] - vp[j]->y, c[2] - vp[j]->z, c[3] - vp[j]->t);
		if (d < mind)
		{
			mind = d;
			rep = vp[j];
		}
	}
	return rep;
}

CoreParticle*
descendUntilCore(CoreParticle* p)
{
	vector<CoreParticle*> Q(1, p);
	while (Q.empty() == false)
	{
		set<CoreParticle*> S;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* q = Q[i];
			if (q->descendents.empty()) return q;

			for (set<CoreParticle*>::iterator it = q->descendents.begin(); it != q->descendents.end(); it++)
			{
				S.insert(*it);
			}
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
	return NULL;
}

vector<CoreParticle*>
closestDescendents(CoreParticle* p)
{
	vector<CoreParticle*> vn = selectNeighbors(p, p->gen + 1);
	vector<CoreParticle*> vq;
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < vn.size(); ++i)
	{
		CoreParticle* q = vn[i];
		if (q->gen == p->gen && q->id < p->id) continue;

		float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
		if (d < mind)
		{
			mind = d;
			vq.clear();
			vq.push_back(q);
		}
		else if (Abs(d-mind) < 0.001f)
		{
			vq.push_back(q);
		}
	}
	return vq;
}

vector<CoreParticle*>
closestInSameGeneration(CoreParticle* p)
{
	vector<CoreParticle*> vq;
	vector<CoreParticle*> vn = selectNeighbors(p, p->gen);
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < vn.size(); ++i)
	{
		CoreParticle* q = vn[i];
		//if (q->gen == p->gen && q->id < p->id) continue;
		for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
		{
			CoreParticle* r = *it;
			float d = length(2 * p->x - r->x - q->x, 2 * p->y - r->y - q->y, 2 * p->z - r->z - q->z, 2 * p->t - r->t - q->t);
			if (d < mind)
			{
				mind = d;
				vq.clear();
				vq.push_back(q);
			}
			else if (Abs(d - mind) < 0.001f)
			{
				vq.push_back(q);
			}
		}
	}
	for (int i = vq.size() - 1; i >= 0; i--)
	{
		if (vq[i]->id < p->id)
		{
			vq.erase(vq.begin() + i);
		}
	}
	return vq;
}

/*
Link particles without descendents to its neighbors with the same generation that is strong medial particle.
This is to link particles inside isolated local maximum components.
*/
void
propagateDescendency2(vector<CoreParticle*>& P, int ndim)
{
	vector<CoreParticle*> U; //particles that can links to decendency
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->descendents.empty())
		{
			U.push_back(P[i]);
		}
	}
	vector<vector<CoreParticle*>> C = clusterParticles(U);
	vector<CoreParticle*> R;
	for (int i = 0; i < C.size(); ++i)
	{
		CoreParticle* rep = closestParticle(C[i], centroid(C[i]));
		if (rep)
		{
			rep->descendents.insert(rep); //temporary add itself as a descendent
			R.push_back(rep);
		}
	}
	vector<CoreParticle*> Q = R;
	while (Q.empty() == false)
	{
		set<CoreParticle*> S; //particles that are to be lined to decendency
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			for (int j = 0; j < p->neighbors8.size(); ++j)
			{
				CoreParticle* q = p->neighbors8[j];
				if (p->gen == q->gen)
				{
					if (q->descendents.size() == 0)
					{
						S.insert(q);
					}
				}
			}
		}
		//for each particle in S, find the best one according to representative neighbor
		for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			CoreParticle* p = *it;
			vector<CoreParticle*> vp;
			vector<CoreParticle*> vn = selectNeighbors(p, p->gen);
			float mind = std::numeric_limits<float>::infinity();
			for (int j = 0; j < vn.size(); ++j)
			{
				CoreParticle* q = vn[j];
				if (q->descendents.empty() == false)
				{
					float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
					if (d < mind)
					{
						vp.clear();
						vp.push_back(q);
						mind = d;
					}
					else if (Abs(d - mind) < 0.001f)
					{
						vp.push_back(q);
					}
				}
			}
			for (int j = 0; j < vp.size(); ++j)
			{
				CoreParticle* q = vp[j];
				p->descendents.insert(q);
				q->ascendents.insert(p);
			}
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
	for (int i = 0; i < R.size(); ++i)
	{
		R[i]->descendents.clear();
	}
}

/*
Link particles without descendents to its neighbors with the same generation and with a decendent.
*/
void
propagateDescendency(vector<CoreParticle*>& P)
{
	vector<CoreParticle*> Q; //particles that can links to decendency
	for (int i = 0; i < P.size(); ++i)
	{
		if (P[i]->descendents.size()>0)
		{
			Q.push_back(P[i]);
		}
	}
	while (Q.empty() == false)
	{
		set<CoreParticle*> S; //particles that are to be linked to decendency
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			for (int j = 0; j < p->neighbors8.size(); ++j)
			{
				CoreParticle* q = p->neighbors8[j];
				if (p->gen == q->gen)
				{
					if (q->descendents.size() == 0)
					{
						S.insert(q);
					}
				}
			}
		}
		//for each particle in S, find the best one according to representative neighbor
		for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			CoreParticle* p = *it;
			vector<CoreParticle*> vp;
			vector<CoreParticle*> vn = selectNeighbors(p, p->gen);
			float mind = std::numeric_limits<float>::infinity();
			for (int j = 0; j < vn.size(); ++j)
			{
				CoreParticle* q = vn[j];
				if (q->descendents.empty() == false)
				{
					float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
					if (d < mind)
					{
						vp.clear();
						vp.push_back(q);
						mind = d;
					}
					else if (Abs(d - mind) < 0.001f)
					{
						vp.push_back(q);
					}
				}
			}
			for (int j = 0; j < vp.size(); ++j)
			{
				CoreParticle* q = vp[j];
				p->descendents.insert(q);
				q->ascendents.insert(p);
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
	vector<CoreParticle*> Q;
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		p->value = 0;
		if (surfaceParticle(p, ndim))
		{
			Q.push_back(p);
		}
	}
	printf("There are %d surface voxels.\n", Q.size());
	int gen = 1;
	while (Q.empty() == false)
	{
		//vector<CoreParticle*> core;
		for (int i = 0; i < Q.size(); ++i)
		{
			SetVoxel(S, Q[i], gen, ndim, dims);
			Q[i]->gen = gen;
			Q[i]->value = 2;
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
					q->value = 1;
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());

		gen++;
	}
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		if (_debug_id == p->id)
		{
			_debug_id += 0;
		}
		vector<CoreParticle*> vq = closestDescendents(p);
		for (int j = 0; j < vq.size(); ++j)
		{
			p->descendents.insert(vq[j]);
			vq[j]->ascendents.insert(p);
		}
	}
	map<CoreParticle*, vector<CoreParticle*>> msp;
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		if (p->descendents.empty())
		{
			vector<CoreParticle*> vq = closestInSameGeneration(p);
			if (vq.empty() == false)
			{
				msp[p] = vq;
			}
		}
	}
	for (map<CoreParticle*, vector<CoreParticle*>>::iterator it = msp.begin(); it != msp.end(); ++it)
	{
		CoreParticle* p = (it)->first;
		vector<CoreParticle*> vq = it->second;
		for (int j = 0; j < vq.size(); ++j)
		{
			p->descendents.insert(vq[j]);
			vq[j]->ascendents.insert(p);
		}
	}

	//propagateDescendency(particles);
	//propagateDescendency2(particles, ndim);

	//remove non-essential link
	/*while (prune)
	{
		bool bChanged = false;
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			if (p->ascendents.size()>1)
			{
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->descendents.size() > 1)
					{
						p->ascendents.erase(find(p->ascendents.begin(), p->ascendents.end(), q));
						q->descendents.erase(find(q->descendents.begin(), q->descendents.end(), p));
						bChanged = true;
						break;
					}
				}
			}
		}
		if (!bChanged) break;
	}*/

	return particles;
}

vector<Vertex<CoreParticle*>*>
makeGraphStructure(vector<CoreParticle*>& mp, vector<int>& S, 
				   vector<CoreParticle*>& selected, //these points will be added to the skeleton.
				   int ndim, const int* dims)
{
	vector<CoreParticle*> core;
	vector<CoreParticle*> particles;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			mp[i]->label = 0;
			mp[i]->value = 0;
			particles.push_back(mp[i]);
		}
	}
	labelStrongCoreParticles(particles, ndim);

	//make sure selected points are core.
	selected[0]->core = 3;
	selected[1]->core = 4;

	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i]->core > 0)
		{
			core.push_back(particles[i]);
			particles[i]->src = particles[i];
		}
		else {
			particles[i]->src = NULL;
		}
	}

	propagateDescendency(particles);
	propagateDescendency2(particles, ndim);

	set<int> sgen;
	for (int i = 0; i < core.size(); ++i)
	{
		sgen.insert(core[i]->gen);
		core[i]->label = i + 1;
	}
	vector<int> vgen;
	vgen.insert(vgen.begin(), sgen.begin(), sgen.end());
	sort(vgen.begin(), vgen.end());
	for (int ig = 0; ig < vgen.size(); ++ig)
	{
		int gval = vgen[ig];
		vector<CoreParticle*> Q;
		for (int i = 0; i < core.size(); ++i)
		{
			if (core[i]->gen == gval)
			{
				Q.push_back(core[i]);
			}
		}
		while (Q.empty() == false)
		{
			set<CoreParticle*> Q2;
			for (int i = 0; i < Q.size(); ++i)
			{
				CoreParticle* p = Q[i];
				SetVoxel(S, p, p->label, ndim, dims);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
					//for (vector<CoreParticle*>::iterator it = p->neighbors.begin(); it != p->neighbors.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->label <= 0)
					{
						q->label = p->label;
						q->value = p->value + 1;
						q->src = p->src;
						Q2.insert(q);
					}
				}
			}
			Q.clear();
			Q.insert(Q.end(), Q2.begin(), Q2.end());
		}
	}
	//now find labels that are adjacent to each other
	GraphFactory<CoreParticle*>& factory = GraphFactory<CoreParticle*>::GetInstance();
	vector<Vertex<CoreParticle*>*> vertices;
	for (int i = 0; i < core.size(); ++i)
	{
		Vertex<CoreParticle*>* u = factory.makeVertex(core[i]);
		vertices.push_back(u);
		core[i]->vertex = u;
	}
	vector<Edge<CoreParticle*>*> edges;
	//map<pair<CoreParticle*, CoreParticle*>, Edge<CoreParticle*>*> emap;
	int edge_count = 0;
	for (int i = 0; i < mp.size(); ++i)
	{
		CoreParticle* p = mp[i];
		if (p)
		{
			for (int j = 0; j < p->neighbors.size(); ++j)
			{
				CoreParticle* q = p->neighbors[j];
				if (p->id == 840 && q->id == 1551)
				{
					j += 0;
				}
				if (p->label > 0 && q->label > 0 && p->label != q->label)
				{
					Vertex<CoreParticle*>* u = vertices[p->label - 1];
					Vertex<CoreParticle*>* v = vertices[q->label - 1];
					Edge<CoreParticle*>* uv = u->findEdge(v);
					Edge<CoreParticle*>* vu = v->findEdge(u);
					float w = weight(p, q);
					if (uv == NULL)
					{
						uv = factory.makeEdge(u, v, w);
						vu = factory.makeEdge(v, u, weight(p, q));
						u->Add(uv);
						v->Add(vu);
						edge_count += 2;
						edges.push_back(uv);
						edges.push_back(vu);
					}
					if (w < uv->w)
					{
						uv->w = w;
						vu->w = w;
					}
					//printf("%d %3.3f -- %d %3.3f ==> %f\n", p->id, p->dval, q->id, q->dval, uv->w);
				}
			}
		}
	}
	//vector<Edge<CoreParticle*>*> mst = Kruskal(vertices);
	//printf("#vertices = %d, #edges = %d\n", vertices.size(), edge_count);

	/*for (int i = 0; i < mst.size(); ++i)
	{
		printf("%d %3.3f -- %d %3.3f ==> %f\n", 
			mst[i]->u->key->id, mst[i]->u->key->dval, mst[i]->v->key->id, mst[i]->v->key->dval, mst[i]->w);
	}*/
	return vertices;
}

void
partition(vector<Edge<CoreParticle*>*>& tree, vector<CoreParticle*>& particles, vector<int>& S, 
		  set<Edge<CoreParticle*>*>& toremove, int ndim, const int* dims)
{
	Vertex<CoreParticle*>* src = NULL;
	set<Vertex<CoreParticle*>*> vertices;
	for (int i = 0; i < tree.size(); ++i) {
		vertices.insert(tree[i]->u);
		vertices.insert(tree[i]->v);
		if (tree[i]->u->key->core == 3)
		{
			src = tree[i]->u;
		}
	}
	assert(src != NULL);

	vector<Node<Vertex<CoreParticle*>*>*> nodes;
	map<Vertex<CoreParticle*>*, Node<Vertex<CoreParticle*>*>*> vnmap;
	for (set<Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
		Node<Vertex<CoreParticle*>*>* n = makeset(*it);
		nodes.push_back(n);
		vnmap[*it] = n;
	}
	for (int i = 0; i < tree.size(); ++i){
		if (toremove.find(tree[i]) == toremove.end())
		{
			merge(vnmap[tree[i]->u], vnmap[tree[i]->v]);
		}
	}
	vector<Node<Vertex<CoreParticle*>*>*> reps = clusters(nodes);
	map<Node<Vertex<CoreParticle*>*>*,int> nimap;
	for (int i = 0; i < reps.size(); ++i)
	{
		nimap[reps[i]] = i + 1;
	}
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
	}
	set<int> sgen;
	for (set<Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it) 
	{
		Node<Vertex<CoreParticle*>*>* n = vnmap[*it];
		Node<Vertex<CoreParticle*>*>* r = findset(n);
		n->key->key->label = nimap[r];
		sgen.insert(n->key->key->gen);
	}
	vector<int> vgen;
	vgen.insert(vgen.begin(), sgen.begin(), sgen.end());
	sort(vgen.begin(), vgen.end());
	for (int ig = 0; ig < vgen.size(); ++ig)
	{
		int gval = vgen[ig];
		vector<CoreParticle*> Q;
		for (set<Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it) 
		{
			if ((*it)->key->gen == gval)
			{
				Q.push_back((*it)->key);
			}
		}
		while (Q.empty() == false)
		{
			set<CoreParticle*> Q2;
			for (int i = 0; i < Q.size(); ++i)
			{
				CoreParticle* p = Q[i];
				SetVoxel(S, p, p->label, ndim, dims);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->label <= 0)
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
	}

	//lastly, re-color the volume from the tree with the source in it
	vector<CoreParticle*> Q;
	Node<Vertex<CoreParticle*>*>* nsrc = findset(vnmap[src]);
	for (int i = 0; i < nodes.size(); ++i)
	{
		if (findset(nodes[i]) == nsrc)
		{
			Q.push_back(nodes[i]->key->key);
		}
	}
	while (Q.empty() == false)
	{
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			SetVoxel(S, p, src->key->label, ndim, dims);
			for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
			{
				CoreParticle* q = *it;
				if (q->label != src->key->label)
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


	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
}

/*
return true if the edge is on the shortest path tree.
*/
bool
onTree(Edge<CoreParticle*>* edge)
{
	return edge->v->pi != NULL && edge->v->pi == edge->u;
}

/*
From a given edge, collect edges that follows it on a (directed) tree.
*/
vector<Edge<CoreParticle*>*>
collectBranch(Edge<CoreParticle*>* root_edge)
{
	vector<Edge<CoreParticle*>*> branch;
	vector<Edge<CoreParticle*>*> Q(1, root_edge);
	while (Q.empty() == false)
	{
		vector<Edge<CoreParticle*>*> Q2;
		for (int j = 0; j < Q.size(); ++j)
		{
			branch.push_back(Q[j]);
			for (int k = 0; k < Q[j]->v->aList.size(); ++k)
			{
				if (onTree(Q[j]->v->aList[k]))
				{
					Q2.push_back(Q[j]->v->aList[k]);
				}
			}
		}
		Q = Q2;
	}
	return branch;
}

/*
Find how a branch is deviated from the main branch (where the src-sink path is included).
The calculation takes every leaf node in the branch and find the closest vertex in the main branch.
The function returns the maximum of the distance among the leaf nodes.
*/
float separationFromMainBranch(vector<Edge<CoreParticle*>*>& branch, vector<Edge<CoreParticle*>*>& main_branch)
{
	float maxd = 0;
	for (int i = 0; i < branch.size(); ++i)
	{
		Vertex<CoreParticle*>* v = branch[i]->v;
		//check if v is a leaf
		bool bleaf = true;
		for (int j = 0; j < v->aList.size(); ++j)
		{
			if (v->aList[j]->v->pi == v)
			{
				bleaf = false;
				break;
			}
		}
		if (bleaf)
		{
			CoreParticle* p = branch[i]->v->key;
			for (int j = 0; j < main_branch.size(); ++j)
			{
				CoreParticle* q = main_branch[j]->v->key;
				float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
				maxd = Max(d, maxd);
			}
		}
	}
	return maxd;
}

/*
Find how a branch is deviated from the main branch (where the src-sink path is included).
The calculation takes every node in the branch and find the distance to the sink node.
The function returns the minimum of the distance among the nodes.
*/
float separationFromMainBranch2(vector<Edge<CoreParticle*>*>& branch, Vertex<CoreParticle*>* sink)
{
	CoreParticle* q = sink->key;
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < branch.size(); ++i)
	{
		Vertex<CoreParticle*>* v = branch[i]->v;
		CoreParticle* p = v->key;
		float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
		mind = Min(d, mind);
	}
	return mind;
}

/*
Calculate the cosine of the angle formed by a-b-c.
*/
float 
cosine(CoreParticle* a, CoreParticle* b, CoreParticle* c)
{
	float dp = dotProduct(a->x - b->x, a->y - b->y, a->z - b->z, a->t - b->t, c->x - b->x, c->y - b->y, c->z - b->z, c->t - b->t);
	float ab = length(a->x - b->x, a->y - b->y, a->z - b->z, a->t - b->t);
	float cb = length(c->x - b->x, c->y - b->y, c->z - b->z, c->t - b->t);
	return dp / (ab*cb);
}

/*
Find how a branch is deviated from the main branch (where the src-sink path is included).
The calculation takes every leaf node in the branch and calcuates cosine of the angle of leaf-src-sink.
The function returns the median of the distance among the leaf nodes.
*/
float separationFromMainBranch3(vector<Edge<CoreParticle*>*>& branch, Vertex<CoreParticle*>* src, Vertex<CoreParticle*>* sink)
{
	vector<float> vcos;
	CoreParticle* p = src->key;
	CoreParticle* q = sink->key;
	float len0 = length(q->x - p->x, q->y - p->y, q->z - p->z, q->t - p->t);
	for (int i = 0; i < branch.size(); ++i)
	{
		Vertex<CoreParticle*>* v = branch[i]->v;
		//check if v is a leaf
		bool bleaf = true;
		for (int j = 0; j < v->aList.size(); ++j)
		{
			if (v->aList[j]->v->pi == v)
			{
				bleaf = false;
				break;
			}
		}
		if (bleaf)
		{
			CoreParticle* r = branch[i]->v->key;
			float len = length(r->x - p->x, r->y - p->y, r->z - p->z, r->t - p->t);
			float dp = dotProduct(r->x - p->x, r->y - p->y, r->z - p->z, r->t - p->t, q->x - p->x, q->y - p->y, q->z - p->z, q->t - p->t);
			if (len > 0)
			{
				vcos.push_back(dp / (len0 * len));
			}
			//printf("(%d,%d,%d) %3.3f, %3.3f, %3.3f\n", r->x, r->y, r->z, len, dp, dp / (len*len0));
			printf("%d %d %d, %d %d %d => %f\n",
				q->x - p->x, q->y - p->y, q->z - p->z, r->x - p->x, r->y - p->y, r->z - p->z, dp / (len*len0));
		}
	}
	sort(vcos.begin(), vcos.end());
	return vcos[vcos.size()/2];
}

/*
Find how a branch is deviated from the main branch (where the src-sink path is included).
The calculation first computes the distance from each node to the main branch.
A branch is retained if no node in the branch deviated away from the main branch too much.
*/
float separationFromMainBranch4(vector<Edge<CoreParticle*>*>& branch, vector<Edge<CoreParticle*>*>& main_branch)
{
	float maxd = 0;
	for (int i = 0; i < branch.size(); ++i)
	{
		CoreParticle* p = branch[i]->v->key;
		float mind = std::numeric_limits<float>::infinity();
		for (int j = 0; j < main_branch.size(); ++j)
		{
			CoreParticle* q = main_branch[j]->v->key;
			float d = length(p-> x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
			if (d < mind)
			{
				mind = d;
			}
		}
		if (mind > maxd)
		{
			maxd = mind;
		}
	}
	return maxd;
}

/*
Find how a branch is deviated from the main branch (where the src-sink path is included).
This version returns an edge to be cut, instead of the measure of diviation from the main branch.
It assumes the branch is ordered in sequence starting from the source node.
It evaluates each node along the branch how far it is to the main branch, until the measure exceeds the threshold.
The edge led to the node diviating above threshold needs to be cut and is returned from this method.
*/
vector<Edge<CoreParticle*>*>
separationFromMainBranch5(Edge<CoreParticle*>* edge, vector<Edge<CoreParticle*>*>& main_branch, float thres)
{
	vector<Edge<CoreParticle*>*> R;
	CoreParticle* p = edge->v->key;
	float mind = std::numeric_limits<float>::infinity();
	for (int j = 0; j < main_branch.size(); ++j)
	{
		CoreParticle* q = main_branch[j]->v->key;
		float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
		if (d < mind)
		{
			mind = d;
		}
	}
	if (mind > thres)
	{
		R.push_back(edge);
	}
	else
	{
		for (int i = 0; i < edge->v->aList.size(); ++i)
		{
			Edge<CoreParticle*>* edge2 = edge->v->aList[i];
			if (onTree(edge2))
			{
				vector<Edge<CoreParticle*>*> R2 = separationFromMainBranch5(edge2, main_branch, thres);
				R.insert(R.end(), R2.begin(), R2.end());
			}
		}
	}
	return R;
}

/*
From a set of particles (SP), go up the ascendent graph freely and collect all surface particles that are reachable.
*/
set<CoreParticle*>
freeAscent(set<CoreParticle*>& sp, int ndim)
{
	set<CoreParticle*> all;
	CoreParticle* src = *sp.begin();
	for (set<CoreParticle*>::iterator it = sp.begin(); it != sp.end(); ++it)
	{
		(*it)->src = src;
	}
	set<CoreParticle*> res;
	vector<CoreParticle*> Q;
	Q.insert(Q.begin(), sp.begin(), sp.end());
	while (Q.empty() == false)
	{
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			all.insert(Q[i]);
			if (surfaceParticle(Q[i], ndim))
			{
				res.insert(Q[i]);
			}
			else {
				for (set<CoreParticle*>::iterator it = Q[i]->ascendents.begin(); it != Q[i]->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					if (q->src != src)
					{
						Q2.insert(q);
						q->src = src;
					}
				}
			}
		}
		Q.clear();
		Q.insert(Q.begin(), Q2.begin(), Q2.end());
	}
	for (set<CoreParticle*>::iterator it = all.begin(); it != all.end(); ++it)
	{
		(*it)->src = NULL;
	}
	return res;
}

/*
Find how a branch is deviated from the main branch (where the src-sink path is included).
This version returns an edge to be cut, instead of the measure of diviation from the main branch.
It assumes the branch is ordered in sequence starting from the source node.
From each branch edge, collect a sequence of tree edges until another branch.
For the branch, find the surface particles that are reachable, and they find the bounding box intersection between
the branch surface and main-path surface (pre-calcuated in main_bbox).
*/
vector<Edge<CoreParticle*>*>
separationFromMainBranch6(Edge<CoreParticle*>* edge, BoundingBox& main_bbox, float thres, int ndim)
{
	vector<Edge<CoreParticle*>*> R;
	Edge<CoreParticle*>* ed = edge;
	vector<CoreParticle*> path;
	path.push_back(ed->u->key);
	while (true)
	{
		path.push_back(ed->v->key);
		vector<Edge<CoreParticle*>*> next;
		for (int i = 0; i < ed->v->aList.size(); ++i)
		{
			Edge<CoreParticle*>* ed2 = ed->v->aList[i];
			if (onTree(ed2))
			{
				next.push_back(ed2);
			}
		}
		if (next.size() == 1)
		{
			ed = next[0];
		}
		else
		{
			//check for overlap
			set<CoreParticle*> spath;
			spath.insert(path.begin(), path.end());
			set<CoreParticle*> sf = freeAscent(spath, ndim);
			BoundingBox bb(sf);
			float olap = overlapVolume(bb, main_bbox);
			float ratio = olap / bb.volume();
			for (vector<CoreParticle*>::iterator jt = path.begin(); jt != path.end(); ++jt)
			{
				CoreParticle* r = *jt;
				printf("%d(%d,%d,%d), ", r->id, r->x, r->y, r->z);
			}
			printf("\n\t[%d,%d,%d,]==[%d,%d,%d] --- overlap=%f, bb=%f, rate=%f\n", 
				bb.minc.m_X, bb.minc.m_Y, bb.minc.m_Z, bb.maxc.m_X, bb.maxc.m_Y, bb.maxc.m_Z,
				olap, bb.volume(), ratio);
			if (ratio < thres)
			{
				R.push_back(edge);
			}
			else
			{
				for (int i = 0; i < ed->v->aList.size(); ++i)
				{
					Edge<CoreParticle*>* ed2 = ed->v->aList[i];
					if (onTree(ed2))
					{
						vector<Edge<CoreParticle*>*> R2 = separationFromMainBranch6(ed2, main_bbox, thres, ndim);
						R.insert(R.end(), R2.begin(), R2.end());
					}
				}
			}
			break;
		}
	}
	return R;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [S] = makeGraphStructure(L, selection)");
		return;
	}

	//load figure-ground segmentation
	int ndimL;
	const int* dimsL;
	mxClassID classL;
	vector<unsigned char> L;
	LoadData(L, prhs[0], classL, ndimL, &dimsL);

	vector<int> coords; //boundary and internal points of a ventricle to be segmented
	int ndimV;
	const int* dimsV;
	mxClassID classV;
	LoadData(coords, prhs[1], classV, ndimV, &dimsV);
	if (coords.size() < 6)
	{
		mexErrMsgTxt("makeGraphStructure: selection input requires 6 values (for two 3D points).");
		return;
	}
	float cutoff_thres = 100;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		//int value;
		ReadScalar(cutoff_thres, prhs[2], classMode);
	}

	vector<float> vsrc(4, 0.0f); 
	vsrc[0] = coords[0], vsrc[1] = coords[1], vsrc[2] = coords[2];
	vector<float> vsink(4, 0.0f);
	vsink[0] = coords[3], vsink[1] = coords[4], vsink[2] = coords[5];

	int nvoxels = numberOfElements(ndimL, dimsL);

	vector<CoreParticle*> mp = generateParticleMap(L, ndimL, dimsL);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndimL, dimsL);
	vector<float> D(L.size(), 0.0f);
	{
		vector<unsigned char> iL(L.size(), 0);
		for (int i = 0; i < L.size(); ++i)
		{
			iL[i] = L[i] ? 0 : 1;
		}
		vector<float> vs(ndimL, 1.0f);
		DistanceTransformEuclidF(D, iL, vs, ndimL, dimsL);
		iL.clear(); //to save memory 
		for (int i = 0; i < D.size(); ++i)
		{
			if (mp[i] != NULL)
			{
				mp[i]->dval = D[i];
			}
		}
	}

	vector<int> S(nvoxels, 0);
	propagateParticles(particles, mp, S,  ndimL, dimsL);

	vector<CoreParticle*> selected_core(2);
	selected_core[0] = closestParticle(particles, vsrc);
	selected_core[1] = closestParticle(particles, vsink);
	//selected_core[0] = descendUntilCore(selected_core[0]);
	//selected_core[1] = descendUntilCore(selected_core[1]);

	vector<int> S2(nvoxels, 0);
	vector<Vertex<CoreParticle*>*> vertices = makeGraphStructure(mp, S2, selected_core, ndimL, dimsL);
	Vertex<CoreParticle*>* src = selected_core[0]->vertex;
	Vertex<CoreParticle*>* sink = selected_core[1]->vertex;

	Dijkstra(vertices, src);
	vector<Edge<CoreParticle*>*> tree;
	for (int i = 0; i < vertices.size(); ++i)
	{
		if (vertices[i]->pi != NULL)
		{
			Edge<CoreParticle*>* ed = vertices[i]->pi->findEdge(vertices[i]);
			assert(ed != NULL);
			tree.push_back(ed);
			ed->type = (EdgeType) 1;
		}
	}

	//clear src pointer
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->src = NULL;
	}

	Vertex<CoreParticle*>* u = sink;
	vector<Edge<CoreParticle*>*> path;
	set<CoreParticle*> spath;
	while (u != src && u != NULL)
	{
		spath.insert(u->key);
		if (u->pi != NULL)
		{
			Edge<CoreParticle*>* ed = u->pi->findEdge(u);
			path.insert(path.begin(), ed);
			ed->type = (EdgeType)2;
		}
		u = u->pi;
	}
	spath.insert(src->key);
	BoundingBox main_bbox(freeAscent(spath, ndimL));
	printf("MainBox: [%d,%d,%d]==[%d,%d,%d]\n",
		main_bbox.minc.m_X, main_bbox.minc.m_Y, main_bbox.minc.m_Z,
		main_bbox.maxc.m_X, main_bbox.maxc.m_Y, main_bbox.maxc.m_Z);
	vector<Edge<CoreParticle*>*> main_branch;
	Edge<CoreParticle*>* retain;
	for (int i = 0; i < src->aList.size(); ++i)
	{
		Edge<CoreParticle*>* ed = src->aList[i];
		if (onTree(ed))
		{
			if (find(path.begin(), path.end(), src->aList[i]) != path.end())
			{
				main_branch = collectBranch(src->aList[i]);
				retain = src->aList[i];
				break;
			}
		}
	}

	/*{
		//for debugging
		Vertex<CoreParticle*>* u = NULL;
		for (int i = 0; i < vertices.size(); ++i)
		{
			if (vertices[i]->key->id == 14275)
			{
				u = vertices[i];
				break;
			}
		}
		printf("%x\n", src->pi);
		src->pi = NULL;
		while (u->pi != NULL)
		{
			set<CoreParticle*> sp;
			sp.insert(u->key);
			sp.insert(u->pi->key);
			set<CoreParticle*> sp2 = freeAscent(sp, ndimL);
			BoundingBox bb(sp2);
			float olap = overlapVolume(bb, main_bbox);
			printf("DBG:: %d (%d, %d, %d) - %d (%d, %d, %d), %f, %f\n", 
				u->key->id, u->key->x, u->key->y, u->key->z, 
				u->pi->key->id, u->pi->key->x, u->pi->key->y, u->pi->key->z,
				olap, bb.volume()); 
			u = u->pi;
		}
	}*/

	/*//for testing
	for (int i = 0; i < path.size(); ++i)
	{
		//vector<Edge<CoreParticle*>*> cut = separationFromMainBranch6(path[i], main_bbox, cutoff_thres, ndimL);
		set<CoreParticle*> sp;
		sp.insert(path[i]->u->key);
		sp.insert(path[i]->v->key);
		set<CoreParticle*> spa = freeAscent(sp, ndimL);
		BoundingBox bb(spa);
		float ol = overlapVolume(bb, main_bbox);
		printf("path (%d,%d,%d) - (%d,%d,%d), %f, %f\n", path[i]->u->key->x, path[i]->u->key->y, path[i]->u->key->z,
			path[i]->v->key->x, path[i]->v->key->y, path[i]->v->key->z,
			ol, bb.volume());
	}*/

	set<Edge<CoreParticle*>*> toremove;
	printf("src = (%d,%d,%d), sink=(%d,%d,%d)\n", src->key->x, src->key->y, src->key->z, sink->key->x, sink->key->y, sink->key->z);
	vector<Edge<CoreParticle*>*> edge_list = src->aList;
	//edge_list.insert(edge_list.end(), retain->v->aList.begin(), retain->v->aList.end());
	for (int i = 0; i < edge_list.size(); ++i)
	{
		Edge<CoreParticle*>* ed = edge_list[i];
		if (onTree(ed))// && ed != retain)
		{
			vector<Edge<CoreParticle*>*> branch = collectBranch(ed);
			//float d = separationFromMainBranch(branch, main_branch);
			//float d = separationFromMainBranch2(branch, sink);
			//float d = separationFromMainBranch3(branch, sink, src);
			//float d = -cosine(ed->v->key, ed->u->key, path[path.size()-1]->u->key);
			/*float d = separationFromMainBranch4(branch, main_branch);
			printf("branch %d - %d: d=%f\n", ed->u->key->id, ed->v->key->id, d);
			if (d > cutoff_thres)
			{
				toremove.insert(ed);
				ed->type = (EdgeType)0;
			}*/
			//vector<Edge<CoreParticle*>*> cut = separationFromMainBranch5(ed, path, cutoff_thres);
			if (ed == retain)
			{
				i += 0;
			}
			vector<Edge<CoreParticle*>*> cut = separationFromMainBranch6(ed, main_bbox, cutoff_thres, ndimL);
			if (ed == retain)
			{
				printf("main branch (%d, %d):: %d to cut\n", path.size(), branch.size(), cut.size());
			}
			for (int k = 0; k < cut.size(); ++k)
			{
				printf("Edge (%d - %d) removed.\n", cut[k]->u->key->id, cut[k]->v->key->id);
				toremove.insert(cut[k]);
				cut[k]->type = (EdgeType)0;
			}
		}
	}

	vector<int> S3(nvoxels, 0);
	partition(tree, particles, S3, toremove, ndimL, dimsL);
	//colorAscendents(particles, S3, toremove, ndimL, dimsL);

	printf("Label at sink = %d\n", sink->key->label);
	if (nlhs >= 1)
	{
		int dims[] = { particles.size(), 12 };
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
			SetData2(F, i, 8, dims[0], dims[1], p->core);
			SetData2(F, i, 9, dims[0], dims[1], inflectionParticle(p, ndimL) ? (int)p->descendents.size() : 0);
			SetData2(F, i, 10, dims[0], dims[1], (int)p->descendents.size());
			SetData2(F, i, 11, dims[0], dims[1], (int)p->ascendents.size());
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { tree.size(), 4 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < tree.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], tree[i]->u->key->id);
			SetData2(F, i, 1, dims[0], dims[1], tree[i]->v->key->id);
			SetData2(F, i, 2, dims[0], dims[1], (int)(100 * tree[i]->w));
			SetData2(F, i, 3, dims[0], dims[1], (int)(tree[i]->type));
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
			//if (surfaceParticle(p, ndimL) == false) continue; //keep it only for surface particles
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
	if (nlhs >= 6)
	{
		plhs[5] = StoreData(S3, mxINT32_CLASS, ndimL, dimsL);
	}
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}

