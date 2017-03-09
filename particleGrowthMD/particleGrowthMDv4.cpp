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
		value2 = 0;
		pi = NULL;
		visited = false;
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
	float value2; //generic value
	bool visited;
	//set<CoreParticle*> ascendents;
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
						//p->ascendents.insert(q);
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
				//p2->ascendents.insert(p);
			}
		}
		propagateDescendency(Q);

		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
		gen++;
	}
	return particles;
}


float measure(vector<CoreParticle*>& particles, vector<CoreParticle*>& mp, 
			  int label1, int label2, int ndim, const int* dims)
{
	vector<float> sum(4, 0.0f);
	for (int i = 0; i < particles.size(); ++i)
	{
		sum[0] += particles[i]->x;
		sum[1] += particles[i]->y;
		sum[2] += particles[i]->z;
		sum[3] += particles[i]->t;
	}
	//get the centroid
	vector<float> mv(4);
	for (int i = 0; i < sum.size(); ++i)
	{
		mv[i] = sum[i] / particles.size();
	}
	//get a particle closest to the centroid
	CoreParticle* src = NULL;
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < particles.size(); ++i)
	{
		float dx = mv[0] - particles[i]->x;
		float dy = mv[1] - particles[i]->y;
		float dz = mv[2] - particles[i]->z;
		float dt = mv[3] - particles[i]->t;
		float dd = sqrt(dx*dx + dy*dy + dz*dz + dt*dt);
		if (dd < mind)
		{
			mind = dd;
			src = particles[i];
		}
		particles[i]->visited = false;
		particles[i]->value2 = 0;
	}
	//initialize particles
	vector<CoreParticle*> Q(1, src);
	src->value2 = 0;
	src->visited = true;
	float maxd = 0.0f;
	float maxval = 0.0f;
	while (Q.empty() == false)
	{
		set<CoreParticle*> S;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			for (int j = 0; j < p->neighbors.size(); ++j)
			{
				CoreParticle* q = p->neighbors[j];
				if ((q->label == label1 || q->label == label2) && q->visited == false)
				{
					S.insert(q);
					q->value2 = p->value2 + 1;
					q->visited = true;
					float df = Abs(src->x - q->x) + Abs(src->y - p->y) + Abs(src->z - q->z) + Abs(src->t - q->t);
					float disc = Abs(df - q->value2);
					if (disc > maxd)
					{
						maxd = disc;
					}
					maxval = Max(maxval, q->value2);
				}
			}
		}
		Q.clear();
		Q.insert(Q.end(), S.begin(), S.end());
	}
	printf("measure: %d - %d w/ %d => %f %f\n", label1, label2, particles.size(), maxval, maxd);
	return maxval >= 0 ? maxd / maxval : 0.0f;
}

/*
void
colorParticles(vector<CoreParticle*>& particles, float thres, int ndim, const int* dims)
{
	vector<CoreParticle*> core;
	vector<CoreParticle*> mp(numberOfElements(ndim, dims), NULL);
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = -1;
		particles[i]->value = std::numeric_limits<float>::infinity();
		SetVoxel(mp, particles[i], particles[i], ndim, dims);
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
		}
		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());
	}
	//find adjacent regions
	vector<vector<bool>> adj;
	for (int i = 0; i < groups.size(); ++i)
	{
		adj.push_back(vector<bool>(groups.size(), false));
	}
	vector<vector<CoreParticle*>> colored(groups.size());
	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		colored[p->label].push_back(p);
		for (int j = 0; j < p->neighbors.size(); ++j)
		{
			CoreParticle* q = p->neighbors[j];
			if (p->label != q->label)
			{
				adj[p->label][q->label] = true;
				adj[q->label][p->label] = true;
			}
		}
	}
	vector<float> isolatedMeasure(groups.size());
	for (int i = 0; i < colored.size(); ++i)
	{
		isolatedMeasure[i] = measure(colored[i], mp, colored[i][0]->label, colored[i][0]->label, ndim, dims);
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
				vector<CoreParticle*> P = colored[i];
				P.insert(P.end(), colored[j].begin(), colored[j].end());
				float mes = measure(P, mp, colored[i][0]->label, colored[j][0]->label, ndim, dims);
				if (mes < (isolatedMeasure[i] + isolatedMeasure[j]) / 2.0f)
				{
					//merge(nodes[i], nodes[j]); ///TK - disable clustering for testing
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
}*/
	
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
	//colorParticles(particles, costCutoff, ndimL, dimsL);

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

