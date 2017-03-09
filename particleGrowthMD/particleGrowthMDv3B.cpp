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
	CoreParticle(int x0 = 0, int y0 = 0, int z0 = 0, int t0 = 0, int gen0 = 0, int lb = 0, int id0 = 0)
	{
		x = x0;
		y = y0;
		gen = gen0;
		label = lb;
		id = id0;
		z = z0;
		t = t0;
		dval = 0;
		value = 0;
		pi = NULL;
	}
	int x;
	int y;
	int z;
	int t;
	int label;
	int gen; //generation
	int id;
	float dval; //distance value
	float value; //generic value
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

vector<CoreParticle*>
setupParticleMap(vector<unsigned char>& L,
vector<float >& D,
int ndim,
const int* dims)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> particles;
	vector<CoreParticle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	vector<CoreParticle*> Q;
	for (int i = 0; i < L.size(); ++i)
	{
		if (L[i])
		{
			vector<int> sub = Ind2Sub(i, ndim, dims);
			CoreParticle* p = coreParticleNdim(sub, ndim);
			mp[i] = p;
			p->dval = GetVoxel(D, p, 0.0f, ndim, dims);
			particles.push_back(p);
		}
	}
	vector<vector<int>>& nbh = NeighborhoodFactory::getInstance(ndim).neighbor4; //MakeEightNeighborhood(ndim, dims);
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
from each surface particle, move towrad the center and establish ascendent/descendent relation.
*/
vector<CoreParticle*>
propagateParticles(
vector<CoreParticle*>& particles,
vector<int>& S,
int ndim,
const int* dims)
{
	vector<CoreParticle*> Q;
	vector<CoreParticle*> mp(numberOfElements(ndim, dims), NULL);
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		SetVoxel(mp, p, p, ndim, dims);
		if (Abs(p->dval - 1.0f) <= 0.1f)
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
			float dval = 0;
			CoreParticle* p2 = NULL;
			for (int n = 0; n < p->neighbors.size(); ++n)
			{
				CoreParticle* q = p->neighbors[n];
				int sval = GetVoxel(S, q, 0, ndim, dims);
				if (sval == 0)
				{
					Q2.insert(q);
					if (q->dval > dval)
					{
						p2 = q;
						dval = q->dval;
					}
				}
			}
			if (p2)
			{
				p->descendents.insert(p2);
				p2->ascendents.insert(p);
			}
		}
		propagateDescendency(Q);

		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
		gen++;
	}
	return particles;
}

template <class T>
struct SymmetricalMap
{
	SymmetricalMap(int N, T defval)
	{
		size = N;
		values = vector<float>(N*N, defval);
		bset = vector<bool>(N*N, false);
	}
	void set(int x, int y, T val)
	{
		int idx = _index(x, y);
		values[idx] = val;
		bset[idx] = true;
	}
	float get(int x, int y)
	{
		return values[_index(x, y)];
	}
	bool valid(int x, int y)
	{
		return bset[_index(x, y)];
	}
	int _index(int x, int y)
	{
		return  Min(x, y) * size + Max(x, y);
	}
	int size;
	vector<T> values;
	vector<bool> bset;
};

//calculate bonding measure based on the width of the boundary, length between the clusters, and depth of the smaller cluster.
//the clusters are tend to be separated if this measure is small.
//the measure is small if w/d is small and len is large.
float measure(float w, float len1, float len2, float d1, float d2)
{
	return w / (sqrt(d1*d2) * sqrt(len1*len2));
}

void
colorParticles(vector<CoreParticle*>& particles, float thres, int ndim, const int* dims)
{
	vector<CoreParticle*> core;
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = -1;
		particles[i]->value = std::numeric_limits<float>::infinity();
		if (particles[i]->descendents.empty())
		{
			core.push_back(particles[i]);
			particles[i]->value = 0;
		}
	}
	vector<vector<CoreParticle*>> groups = clusterParticles(core);
	for (int i = 0; i < groups.size(); ++i)
	{
		for (int j = 0; j < groups[i].size(); ++j)
		{
			groups[i][j]->label = i;
		}
	}
	SymmetricalMap<float> wpmap(groups.size(), 0.0f); //width between two clusters.
	SymmetricalMap<float> dpmap(groups.size(), 0.0f); //distanace from p to the boundary.
	SymmetricalMap<float> dpmap2(groups.size(), 0.0f); //distanace q to the boundary.
	vector<CoreParticle*> Q = core;
	while (Q.empty() == false)
	{
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
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
			for (int n = 0; n < p->neighbors.size(); ++n)
			{
				CoreParticle* q = p->neighbors[n];
				if (q->label >= 0 && p->label != q->label)
				{
					if (wpmap.get(p->label, q->label) < Max(p->dval, q->dval))
					{
						wpmap.set(p->label, q->label, Max(p->dval, q->dval));
						dpmap.set(p->label, q->label, p->label < q->label ? p->value + 1 : q->value + 1);
						dpmap2.set(p->label, q->label, p->label < q->label ? q->value + 1 : p->value + 1);
					}
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
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
			if (dpmap.valid(i, j))
			{
				float bval = measure(wpmap.get(i, j), dpmap.get(i, j), dpmap2.get(i, j), groups[i][0]->dval, groups[j][0]->dval);
				printf("(%d, %d): w=%f, len=%f, len=%f, d=%f, d=%f, measure=%f\n",
					i + 1, j + 1, wpmap.get(i, j), dpmap.get(i, j), dpmap2.get(i, j), groups[i][0]->dval, groups[j][0]->dval, bval);

				if (bval > thres)
				{
					merge(nodes[i], nodes[j]);
				}
			}
		}
	}
	vector<Node<int>*> reps = clusters(nodes);
	printf("%d clusters reduced to %d.\n", nodes.size(), reps.size());
	map<int, int> remap;
	for (int i = 0; i < nodes.size(); ++i)
	{
		remap[i] = distance(reps.begin(), find(reps.begin(), reps.end(), findset(nodes[i])));
	}
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = remap[particles[i]->label];
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

	//int nclusters = 2;
	float costCutoff = 1.0f;
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

	int nvoxels = numberOfElements(ndimL, dimsL);

	vector<float> D(L.size(), 0.0f);
	vector<unsigned char> iL(L.size(), 0);
	for (int i = 0; i<L.size(); ++i)
	{
		iL[i] = L[i] ? 0 : 1;
	}
	vector<float> vs(ndimL, 1.0f);
	DistanceTransformEuclidF(D, iL, vs, ndimL, dimsL);
	iL.clear(); //to save memory 

	vector<CoreParticle*> particles = setupParticleMap(L, D, ndimL, dimsL);
	vector<int> S(nvoxels, 0);
	particles = propagateParticles(particles, S, ndimL, dimsL);
	colorParticles(particles, costCutoff, ndimL, dimsL);

	vector<CoreParticle*> core;
	vector<int> T(nvoxels, 0);
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		SetVoxel(T, p, p->label+1, ndimL, dimsL);
		if (p->descendents.empty())
		{
			core.push_back(p);
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
		int dims[] = { core.size(), 7 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < core.size(); ++i)
		{
			CoreParticle* p = core[i];
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
		plhs[3] = StoreData(T, mxINT32_CLASS, ndimL, dimsL);
	}
	if (nlhs >= 5)
	{
		plhs[4] = StoreData(D, mxSINGLE_CLASS, ndimL, dimsL);
	}
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}

