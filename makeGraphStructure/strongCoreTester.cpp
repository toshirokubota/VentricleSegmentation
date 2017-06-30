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
	float len = length(p->src->x - q->src->x, p->src->y - q->src->y, p->src->z - q->src->z, p->src->t - q->src->t);
	return (len + sqrt(p->src->dval * q->src->dval) + 1.0) / sqrt(p->dval * q->dval);
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
		if (p->ascendents.size() >= (ndim - 1) * 2)
		{
			set<CoreParticle*> S;
			vector<CoreParticle*> Q;
			Q.insert(Q.end(), p->ascendents.begin(), p->ascendents.end());
			while (Q.empty() == false)
			{
				set<CoreParticle*> S2;
				for (int j = 0; j < Q.size(); ++j)
				{
					CoreParticle* q = Q[j];
					if (surfaceParticle(q, ndim))
					{
						S.insert(q);
						q->label = 1;
					}
					for (set<CoreParticle*>::iterator it = q->ascendents.begin(); it != q->ascendents.end(); ++it)
					{
						S2.insert(*it);
					}
				}
				Q.clear();
				Q.insert(Q.end(), S2.begin(), S2.end());
			}
			vector<CoreParticle*> surface;
			//CoreParticle* sp = *(S.begin());
			for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
			{
				CoreParticle* sp = *it;
				for (int i = 0; i < groups.size(); ++i)
				{
					if (find(groups[i].begin(), groups[i].end(), sp) != groups[i].end())
					{
						surface = groups[i];
						break;
					}
				}
			}

			vector<Node<CoreParticle*>*> nodes;
			map<CoreParticle*, Node<CoreParticle*>*> nmap;
			for (int i = 0; i < surface.size(); ++i)
			{
				if (S.find(surface[i]) == S.end())
				{
					Node<CoreParticle*>* n = makeset(surface[i]);
					nodes.push_back(n);
					nmap[surface[i]] = n;
				}
			}

			for (int i = 0; i < nodes.size(); ++i)
			{
				CoreParticle* q = nodes[i]->key;
				for (int j = 0; j < q->neighbors8.size(); ++j)
				{
					CoreParticle* r = q->neighbors8[j];
					if (surfaceParticle(r, ndim) && nmap.find(r) != nmap.end())
					{
						merge(nodes[i], nmap[r]);
					}
				}
			}
			vector<Node<CoreParticle*>*> reps = clusters(nodes);
			int nc = reps.size();
			if (nc > 1)
			{
				p->core += 1;
			}
			printf("%d (%d,%d,%d) with %d clusters: ", p->id, p->x, p->y, p->z, nc);
			vector<int> vclc(nc, 0);
			for (int k = 0; k < nodes.size(); ++k)
			{
				Node<CoreParticle*>* r = findset(nodes[k]);
				int ind = distance(reps.begin(), find(reps.begin(), reps.end(), r));
				vclc[ind]++;
				nodes[k]->key->label = ind + 2;
			}
			for (int k = 0; k < nc; ++k)
			{
				printf("%d ", vclc[k]);
			}
			printf("\n");

			for (int i = 0; i < nodes.size(); ++i)
			{
				delete nodes[i];
			}
		}
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
	for (int i = 0; i < particles.size(); ++i)
	{
		strongMedialParticle(particles[i], groups, ndim);
	}
}

void
labelStrongCoreParticlesTest(vector<CoreParticle*>& particles, CoreParticle* selected, int ndim)
{
	vector<CoreParticle*> surfaces;
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
		if (surfaceParticle(particles[i], ndim))
		{
			surfaces.push_back(particles[i]);
		}
	}
	vector<vector<CoreParticle*>> groups = clusterParticles(surfaces);
	selected->core = -1;
	strongMedialParticle(selected, groups, ndim);
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

vector<Edge<CoreParticle*>*>
makeGraphStructure(vector<CoreParticle*>& mp, vector<int>& S, int ndim, const int* dims)
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
		vertices.push_back(factory.makeVertex(core[i]));
	}
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
	vector<Edge<CoreParticle*>*> mst = Kruskal(vertices);
	printf("#vertices = %d, #edges = %d\n", vertices.size(), edge_count);

	/*for (int i = 0; i < mst.size(); ++i)
	{
		printf("%d %3.3f -- %d %3.3f ==> %f\n", 
			mst[i]->u->key->id, mst[i]->u->key->dval, mst[i]->v->key->id, mst[i]->v->key->dval, mst[i]->w);
	}*/
	return mst;
}

void
partition(vector<Edge<CoreParticle*>*>& tree, vector<CoreParticle*>& particles, vector<int>& S, 
			vector<int>& toremove, int ndim, const int* dims)
{
	for (int i = 0; i < toremove.size(); ++i)
	{
		printf("to-remove id: %d\n", toremove[i]);
	}
	set<int> sids;
	sids.insert(toremove.begin(), toremove.end());

	set<Vertex<CoreParticle*>*> vertices;
	for (int i = 0; i < tree.size(); ++i) {
		vertices.insert(tree[i]->u);
		vertices.insert(tree[i]->v);
	}
	vector<Node<Vertex<CoreParticle*>*>*> nodes;
	map<Vertex<CoreParticle*>*, Node<Vertex<CoreParticle*>*>*> vnmap;
	for (set<Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
		Node<Vertex<CoreParticle*>*>* n = makeset(*it);
		nodes.push_back(n);
		vnmap[*it] = n;
	}
	for (int i = 0; i < tree.size(); ++i){
		if (sids.find(tree[i]->u->key->id) == sids.end() || sids.find(tree[i]->v->key->id) == sids.end())
		{
			merge(vnmap[tree[i]->u], vnmap[tree[i]->v]);
		}
		else
		{
			printf("Edge (%d - %d) removed.\n", tree[i]->u->key->id, tree[i]->v->key->id);
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
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
}


void
colorAscendents(vector<CoreParticle*>& particles, vector<int>& S, vector<int>& srcIds, int ndim, const int* dims)
{
	set<int> sids;
	sids.insert(srcIds.begin(), srcIds.end());

	set<CoreParticle*> src;
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
		if (sids.find(particles[i]->id) != sids.end())
		{
			src.insert(particles[i]);
			particles[i]->label = src.size();
		}
	}
	vector<CoreParticle*> Q;
	Q.insert(Q.end(), src.begin(), src.end());
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
					Q2.insert(q);
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
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

	vector<int> toremove; //id of core particles to be removed
	if (nrhs >= 2)
	{
		int ndimV;
		const int* dimsV;
		mxClassID classV;
		LoadData(toremove, prhs[1], classV, ndimV, &dimsV);
	}
	int selected = -1;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		//int value;
		ReadScalar(selected, prhs[2], classMode);
		//_EightNeighbor = value > 0 ? true : false;
	}
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
	CoreParticle* selectedParticle = NULL;
	propagateParticles(particles, mp, S,  ndimL, dimsL);
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			mp[i]->label = 0;
			mp[i]->value = 0;
			if (mp[i]->id == selected)
			{
				selectedParticle = mp[i]; 
			}
		}
	}
	if (selectedParticle != NULL)
	{
		printf("Examining from %d\n", selectedParticle->id);
		labelStrongCoreParticlesTest(particles, selectedParticle, ndimL);
	}
	vector<int> S2(nvoxels, 0);
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			S2[i] = mp[i]->label;
		}
	}
	vector<int> S3(nvoxels, 0);
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
	if (nlhs >= 2) //not currently used.
	{
		int dims[] = { 1,1 };
		vector<int> F(dims[0] * dims[1]);
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

