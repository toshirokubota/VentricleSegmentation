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
#include <szConvexHull3D.h>

bool _EightNeighbor = false;
bool _Use_New = true;
int _debug_id = 85958; 
int _debug_id2 = 85957;

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
	float dval; //distance value
	bool selected;
	set<CoreParticle*> ascendents;
	set<CoreParticle*> descendents;
	vector<CoreParticle*> neighbors;
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

/*
When two particles (p and q) meet originated from different strong core sduring the growth, and a graph edge will be formed.
This routine provides the weight of the edge.
*/
float weight(CoreParticle* p, CoreParticle* q)
{
	//return sqrt(p->src->dval * q->src->dval) / sqrt(p->dval * q->dval);
	return sqrt(p->src->dval * q->src->dval) / sqrt(p->dval * q->dval) + 1.0 / sqrt(p->dval * q->dval);
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

bool
checkForCollision(vector<vector<int>>& df, int ndim)
{
	if (ndim == 2)
	{
		bool b[8];
		int xoff[] = { -1, 0, 1, -1, 1, -1, 0, 1 };
		int yoff[] = { -1, -1, -1, 0, 0, 1, 1, 1 };
		for (int n = 0; n < 8; ++n)
		{
			b[n] = false;
			for (int i = 0; i < df.size(); ++i)
			{
				if (df[i][0] == xoff[n] && df[i][1] == yoff[n])
				{
					b[n] = true;
					break;
				}
			}
		}
		return (b[0] && b[7]) || (b[1] && b[6]) || (b[2] && b[5]) || (b[3] && b[4]);
	}
	else if (ndim == 3)
	{
		bool b[26];
		int xoff[] = { -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1, -1, 0, 1 };
		int yoff[] = { -1, -1, -1, 0, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 1, 1, 1, -1, -1, -1, 0, 0, 0, 1, 1, 1 };
		int zoff[] = { -1, -1, -1, -1, -1, -1, -1, -1, -1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
		for (int n = 0; n < 26; ++n)
		{
			b[n] = false;
			for (int i = 0; i < df.size(); ++i)
			{
				if (df[i][0] == xoff[n] && df[i][1] == yoff[n] && df[i][2] == zoff[n])
				{
					b[n] = true;
					break;
				}
			}
		}
		return (b[9] && b[11] && b[14] && b[16]) //Z=0 plane, Diagonal, [0 0 1]
			|| (b[10] && b[12] && b[13] && b[15]) //Z=0 plane Horizontal/Vertical, [0
			|| (b[4] && b[10] && b[15] && b[21])  //X=0 plane, Horizontal/Vertical
			|| (b[1] && b[7] && b[18] && b[24])  //X=0 plane, Diagonal
			|| (b[4] && b[12] && b[13] && b[21]) //Y=0 plane, Horizontal/Vertical
			|| (b[3] && b[5] && b[20] && b[22]) //Y=0 plane, Diagonal
			|| (b[0] && b[2] && b[23] && b[25])
			|| (b[1] && b[12] && b[13] && b[24])
			|| (b[6] && b[8] && b[17] && b[19])
			|| (b[7] && b[12] && b[13] && b[18])
			|| (b[0] && b[6] && b[19] && b[25])
			|| (b[3] && b[10] && b[15] && b[22])
			|| (b[2] && b[8] && b[17] && b[23])
			|| (b[5] && b[10] && b[15] && b[20])
			;

	}
}

/*
Ascendents of p are not linearly separable.
*/
bool
strongMedialParticle(CoreParticle* p, int ndim)
{
	if (p->descendents.size() == 0)
	{
		return true;
	}
	if (medialParticle(p, ndim) == false)
	{
		return false;
	}
	if (ndim == 2 || ndim == 3)
	{
		vector<vector<int>> df;
		for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
		{
			CoreParticle* q = *it;
			vector<int> d;
			d.push_back(q->x - p->x);
			d.push_back(q->y - p->y);
			d.push_back(q->z - p->z);
			d.push_back(q->t - p->t);
			df.push_back(d);
		}
		return checkForCollision(df, ndim);
	}
	else
	{
		return false;
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
Based on the L2 distance value, find the neighbor that is most straight.
*/
CoreParticle*
representativeNeighbor(CoreParticle* p, int gen)
{
	float minerr = std::numeric_limits<float>::infinity();
	CoreParticle* best = NULL;
	for (int i = 0; i < p->neighbors.size(); ++i)
	{
		CoreParticle* q = p->neighbors[i];
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
		vector<float> v(4, 0.0f);
		for (int j = 0; j < C[i].size(); ++j)
		{
			v[0] += C[i][j]->x;
			v[1] += C[i][j]->y;
			v[2] += C[i][j]->z;
			v[3] += C[i][j]->t;
		}
		v[0] /= C[i].size(); v[1] /= C[i].size(); v[2] /= C.size(); v[3] /= C.size();
		CoreParticle* rep = NULL;
		float mind = std::numeric_limits<float>::infinity();
		for (int j = 0; j < C[i].size(); ++j)
		{
			float d = (v[0] - C[i][j]->x, v[1] - C[i][j]->y, v[2] - C[i][j]->z, v[3] - C[i][j]->t);
			if (d < mind)
			{
				mind = d;
				rep = C[i][j];
			}
		}
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
				if (q->gen == p->gen && q->descendents.size() > 0)
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
			if (best->id == _debug_id || best->id == _debug_id2)
			{
				printf("C - %d => %d\n", p->id, best->id);
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
Linke\ particles without descendents to its neighbors with the same generation and with a decendent.
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
			if (p->ascendents.empty())
			{
				best = representativeNeighbor(p, p->gen);
			}
			else {
				vector<float> mv(4, 0.f);
				for (set<CoreParticle*>::iterator jt = p->ascendents.begin(); jt != p->ascendents.end(); jt++)
				{
					CoreParticle* r = *jt;
					mv[0] += r->x;
					mv[1] += r->y;
					mv[2] += r->z;
					mv[3] += r->t;
				}
				for (int m = 0; m < 4; ++m)
				{
					mv[m] /= p->ascendents.size();
				}
				float mind = numeric_limits<float>::infinity();
				float predicted[4] = { 2 * p->x - mv[0], 2 * p->y - mv[1], 2 * p->z - mv[2], 2 * p->t - mv[3] };
				for (int j = 0; j < p->neighbors.size(); ++j)
				{
					CoreParticle* q = p->neighbors[j];
					if (q->descendents.size() > 0 && q->gen == p->gen)
					{
						float d = (predicted[0] - q->x)*(predicted[0] - q->x) + (predicted[1] - q->y)*(predicted[1] - q->y) + 
							(predicted[2] - q->z)*(predicted[2] - q->z) + (predicted[3] - q->t)*(predicted[3] - q->t);
						if (d < mind)
						{
							mind = d;
							best = q;
						}
					}
				}
			}
			p->descendents.insert(best);
			best->ascendents.insert(p);
			if (best->id == _debug_id || best->id == _debug_id2)
			{
				printf("C - %d => %d\n", p->id, best->id);
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
		CoreParticle* q = representativeNeighbor(p, p->gen + 1);
		if (q != NULL)
		{
			p->descendents.insert(q);
			q->ascendents.insert(p);
			if (q->id == _debug_id)
			{
				printf("A - %d => %d\n", p->id, q->id);
			}
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
				if (q->id == _debug_id || q->id == _debug_id2)
				{
					printf("B - %d => %d\n", p->id, q->id);
				}
			}
		}
	}
	propagateDescendency(particles);
	propagateDescendency2(particles, ndim);

	return particles;
}

vector<Edge<CoreParticle*>*>
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
				mp[i]->src = mp[i];
			}
			else {
				mp[i]->src = NULL;
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
					Edge<CoreParticle*>* uv = factory.makeEdge(u, v, weight(p, q));
					u->Add(uv);
					Edge<CoreParticle*>* vu = factory.makeEdge(v, u, weight(p, q));
					v->Add(vu);
				}
			}
		}
	}
	
	vector<Edge<CoreParticle*>*> mst = Kruskal(vertices);
	return mst;
}

void
partition(vector<Edge<CoreParticle*>*>& tree, vector<CoreParticle*>& particles, vector<int>& S, 
			int numClusters, int ndim, const int* dims)
{
	vector<pair<float, Edge<CoreParticle*>*>> pairs;
	set<Vertex<CoreParticle*>*> vertices;
	for (int i = 0; i < tree.size(); ++i) {
		pairs.push_back(pair<float, Edge<CoreParticle*>*>(tree[i]->w, tree[i]));
		vertices.insert(tree[i]->u);
		vertices.insert(tree[i]->v);
	}
	sort(pairs.begin(), pairs.end());
	float thres = 0;
	if (numClusters < pairs.size())
	{
		thres = pairs[pairs.size() - numClusters].first;
	}
	vector<Node<Vertex<CoreParticle*>*>*> nodes;
	map<Vertex<CoreParticle*>*, Node<Vertex<CoreParticle*>*>*> vnmap;
	for (set<Vertex<CoreParticle*>*>::iterator it = vertices.begin(); it != vertices.end(); ++it) {
		Node<Vertex<CoreParticle*>*>* n = makeset(*it);
		nodes.push_back(n);
		vnmap[*it] = n;
	}
	for (int i = 0; i < tree.size(); ++i){
		if (tree[i]->w < thres) {
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
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
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

	int nclusters = 2;
	//float cutoff = 3;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(nclusters, prhs[1], classMode);
	}
	if (nrhs >= 3)
	{
		mxClassID classMode;
		int value;
		ReadScalar(value, prhs[2], classMode);
		_EightNeighbor = value > 0 ? true : false;
	}
	if (nrhs >= 4)
	{
		mxClassID classMode;
		int value;
		ReadScalar(value, prhs[3], classMode);
		_Use_New = value > 0 ? true : false;
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
	propagateParticles(particles, mp, S, ndimL, dimsL);
	vector<int> S2(nvoxels, 0);
	vector<Edge<CoreParticle*>*> mst = makeGraphStructure(mp, S2, ndimL, dimsL);
	vector<int> S3(nvoxels, 0);
	partition(mst, particles, S3, nclusters, ndimL, dimsL);

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
			SetData2(F, i, 8, dims[0], dims[1], strongMedialParticle(p, ndimL) ? 1 : 0);
			SetData2(F, i, 9, dims[0], dims[1], inflectionParticle(p, ndimL) ? (int)p->descendents.size() : 0);
			SetData2(F, i, 10, dims[0], dims[1], (int)p->descendents.size());
			SetData2(F, i, 11, dims[0], dims[1], (int)p->ascendents.size());
			if (p->id == 16838 || p->id == 9032)
			{
				bool b = strongMedialParticle(p, ndimL);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					printf("%d(%d,%d,%d) -- %d(%d,%d,%d)\n", p->id, p->x, p->y, p->z, q->id, q->x, q->y, q->z);
				}
				printf("%d %d %d\n", p->x + 1, p->y + 1, p->z + 1);
				for (set<CoreParticle*>::iterator it = p->ascendents.begin(); it != p->ascendents.end(); ++it)
				{
					CoreParticle* q = *it;
					printf("%d %d %d\n", q->x + 1, q->y + 1, q->z + 1);
				}
			}
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { mst.size(), 3 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < mst.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], mst[i]->u->key->id);
			SetData2(F, i, 1, dims[0], dims[1], mst[i]->v->key->id);
			SetData2(F, i, 2, dims[0], dims[1], (int)(100*mst[i]->w));
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

