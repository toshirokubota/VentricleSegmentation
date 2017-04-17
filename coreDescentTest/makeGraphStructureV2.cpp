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
	}
	~NeighborhoodFactory()
	{
	}
	NeighborhoodFactory(NeighborhoodFactory& f){}
	NeighborhoodFactory operator=(NeighborhoodFactory& f){}
};

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

/*
p has ascendents that are not neighbors to each other.
*/
bool
strongMedialParticle(CoreParticle* p, int ndim)
{
	if (medialParticle(p, ndim) == false)
	{
		return false;
	}
	for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
	{
		CoreParticle* q = *it;
		set<CoreParticle*>::iterator jt = it;
		jt++;
		for (; jt != p->ascendents.end(); ++jt)
		{
			CoreParticle* r = *jt;
			if (find(r->neighbors.begin(), r->neighbors.end(), q) == r->neighbors.end()) return true;
		}
	}
	return false;
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
	if (_EightNeighbor)
	{
		nbh = nfactory.neighbor8;
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

/*
Linke particles without descendents to its neighbors with the same generation that is strong medial particle.
*/
void
propagateDescendency2(vector<CoreParticle*>& P, int ndim)
{
	vector<CoreParticle*> Q; //particles that can links to decendency
	for (int i = 0; i < P.size(); ++i)
	{
		if (strongMedialParticle(P[i], ndim))
		{
			Q.push_back(P[i]);
		}
	}
	while (Q.empty() == false)
	{
		set<CoreParticle*> S; //particles that are to be lined to decendency
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
						S.insert(q);
					}
				}
			}
		}
		//for each particle in S, find the best one according to representative neighbor
		for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			CoreParticle* p = *it;
			CoreParticle* best = NULL;
			float mind = numeric_limits<float>::infinity();
			for (int j = 0; j < p->neighbors.size(); ++j)
			{
				CoreParticle* q = p->neighbors[j];
				if (q->gen == p->gen && (q->descendents.size() > 0 || strongMedialParticle(q, ndim)))
				{
					float d = (p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) + (p->z - q->z)*(p->z - q->z) + (p->t - q->t)*(p->t - q->t);
					if (d < mind)
					{
						mind = d;
						best = q;
					}
				}
			}
			p->descendents.insert(best);
			best->ascendents.insert(p);
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
}

/*
Linke particles without descendents to its neighbors with the same generation and with a decendent.
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
		set<CoreParticle*> S; //particles that are to be lined to decendency
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
						S.insert(q);
					}
				}
			}
		}
		//for each particle in S, find the best one according to representative neighbor
		for (set<CoreParticle*>::iterator it = S.begin(); it != S.end(); ++it)
		{
			CoreParticle* p = *it;
			CoreParticle* best = NULL;
			float mind = numeric_limits<float>::infinity();
			for (int j = 0; j < p->neighbors.size(); ++j)
			{
				CoreParticle* q = p->neighbors[j];
				if (q->descendents.size() > 0 && q->gen == p->gen)
				{
					float d = (p->x - q->x)*(p->x - q->x) + (p->y - q->y)*(p->y - q->y) + (p->z - q->z)*(p->z - q->z) + (p->t - q->t)*(p->t - q->t);
					if (d < mind)
					{
						mind = d;
						best = q;
					}
				}
			}
			p->descendents.insert(best);
			best->ascendents.insert(p);
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
}


/*
Based on the surface normal direction, find a neighbor that match the direction best.
GEN can be one larger than p->gen (in this case, searching for a descendent) or
one less than p->gen (in this case, searching for an ancestor.
*/
CoreParticle*
representativeNeighbor(CoreParticle* p, int gen)
{
	float dx = 0, dy = 0, dz = 0, dt = 0;
	int count = 0;
	for (int i = 0; i < p->neighbors.size(); ++i)
	{
		CoreParticle* q = p->neighbors[i];
		if (q->gen == gen)
		{
			dx += (q->x - p->x);
			dy += (q->y - p->y);
			dz += (q->z - p->z);
			dt += (q->t - p->t);
			count++;
		}
	}
	if (count == 0) return NULL;
	else
	{
		float mx = p->x + dx / count;
		float my = p->y + dy / count;
		float mz = p->z + dz / count;
		float mt = p->t + dt / count;
		CoreParticle* best = NULL;
		float mind = std::numeric_limits<float>::infinity();
		for (int i = 0; i < p->neighbors.size(); ++i)
		{
			CoreParticle* q = p->neighbors[i];
			if (q->gen ==  gen)
			{
				float d = (q->x - mx)*(q->x - mx) + (q->y - my)*(q->y - my) + (q->z - mz)*(q->z - mz) + (q->t - mt)*(q->t - mt);
				if (d < mind)
				{
					mind = d;
					best = q;
				}
			}
		}
		return best;
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
		CoreParticle* q = representativeNeighbor(p, p->gen + 1);
		if (q != NULL)
		{
			p->descendents.insert(q);
			q->ascendents.insert(p);
		}
	}
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		if (p->ascendents.empty())
		{
			CoreParticle* q = representativeNeighbor(p, p->gen - 1);
			if (q != NULL)
			{
				q->descendents.insert(p);
				p->ascendents.insert(q);
			}
		}
	}
	propagateDescendency(particles);
	propagateDescendency2(particles, ndim);

	return particles;
}

set<pair<CoreParticle*,CoreParticle*>>
makeGraphStructure(vector<CoreParticle*>& mp, vector<int>& S, int ndim, const int* dims)
{
	vector<CoreParticle*> core;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i])
		{
			mp[i]->label = 0;
			mp[i]->value = 0;
			if (strongMedialParticle(mp[i], ndim))
			{
				core.push_back(mp[i]);
			}
		}
	}
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
					Edge<CoreParticle*>* uv = factory.makeEdge(u, v, p->value + q->value);
					u->Add(uv);
					Edge<CoreParticle*>* vu = factory.makeEdge(v, u, p->value + q->value);
					v->Add(vu);
				}
			}
		}
	}
	
	vector<Edge<CoreParticle*>*> mst = Kruskal(vertices);
	set<pair<CoreParticle*, CoreParticle*>> adj;
	for (int i = 0; i < mst.size(); ++i)
	{
		CoreParticle* p = mst[i]->u->key;
		CoreParticle* q = mst[i]->v->key;
		adj.insert(pair<CoreParticle*, CoreParticle*>(p, q));
	}
	return adj;
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
	set<pair<CoreParticle*, CoreParticle*>> adj = makeGraphStructure(mp, S2, ndimL, dimsL);

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
			SetData2(F, i, 8, dims[0], dims[1], strongMedialParticle(p, ndimL) ? 1 : 0);
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
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}

