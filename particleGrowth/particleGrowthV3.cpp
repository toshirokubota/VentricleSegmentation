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

struct Particle 
{
	Particle(int x0 = 0, int y0 = 0, int gen0=0, int lb=0, int id0 = 0)
	{
		x = x0;
		y = y0;
		gen = gen0;
		label = lb;
		id = id0;
		corelink = NULL;
	}
	int x;
	int y;
	int label;
	int gen; //generation
	int id;
	set<Particle*> ascendents;
	set<Particle*> descendents;
	Particle* corelink; //used to link cores

};

struct _particle_equal{
	_particle_equal(const Particle* p0) : p(p0) { }
	bool operator()(Particle* q) const { return  p->x == q->x && p->y == q->y; }
private:
	const Particle* p;
};

struct ParticleFactory
{
public:
	static ParticleFactory& getInstance()
	{
		static ParticleFactory instance;
		return instance;
	}
	Particle* makeParticle(int x, int y, int gen,  int lb=0)
	{
		Particle* particle = new Particle(x, y, gen, lb, _id++);
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
	vector<Particle*> particles;
private:
	int _id;
	ParticleFactory()
	{
		_id = 0;
	}
	~ParticleFactory()
	{
		clean();
	}
	ParticleFactory(ParticleFactory& f){}
	ParticleFactory operator=(ParticleFactory& f){}
};

const int XOffsetHere[] = { -1, 0, 1, -1, 1, -1, 0, 1 };
const int YOffsetHere[] = { -1, -1, -1, 0, 0, 1, 1, 1 };
const int NumNeighborsHere = sizeof(XOffsetHere) / sizeof(XOffsetHere[0]);
const int XOffsetFour[] = { 0, -1, 1, 0 };
const int YOffsetFour[] = { -1, 0, 0, 1 };
const int NumNeighborsFour = sizeof(XOffsetFour) / sizeof(XOffsetFour[0]);
const int XOffsetDiagonal[] = { -1, 1, -1, 1 };
const int YOffsetDiagonal[] = { -1, -1, 1, 1 };
const int NumNeighborsDiagonal = sizeof(XOffsetFour) / sizeof(XOffsetFour[0]);

template<class T>
bool
SetVoxel(vector<T>& A,
const Particle* p,
T value,
const int* dims)
{
	return SetData2(A, p->x, p->y, dims[0], dims[1], value);
}

template<class T>
T
GetVoxel(const vector<T>& A,
const Particle* p,
T defaultValue,
const int* dims)
{
	return GetData2(A, p->x, p->y, dims[0], dims[1], defaultValue);
}

template<class T>
bool
SetVoxel(vector<T>& A,
const CParticle4D& p,
T value,
const int* dims)
{
	return SetData2(A, p.m_X, p.m_Y, dims[0], dims[1], value);
}

template<class T>
T
GetVoxel(const vector<T>& A,
const CParticle4D& p,
T defaultValue,
const int* dims)
{
	return GetData2(A, p.m_X, p.m_Y, dims[0], dims[1], defaultValue);
}

/*
onSurface3 with a general neighborhood
*/
template<class T>
bool
onSurface2(const vector<T>& A, int x, int y, const int* dims) 
{
	if (GetData2(A, x, y, dims[0], dims[1], (T)0))
	{
		for(int m=0; m<NumNeighborsHere; ++m)
		{
			int x2 = x + XOffsetHere[m];
			int y2 = y + YOffsetHere[m];
			if (GetData2(A, x2, y2, dims[0], dims[1], (T)1) == 0)
			{
				return true;
			}
		}
	}
	return false;
}

vector<Particle*>
ExtactSurface(const vector<unsigned char>& L, 
			  const int* dims)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<Particle*> front;
	for (int j = 0; j < dims[1]; ++j)
	{
		for (int k = 0; k < dims[0]; ++k)
		{
			if (onSurface2(L, k, j, dims))
			{
				front.push_back(factory.makeParticle(k, j, 0));
			}
		}
	}
	return front;
}

float 
Length2(float x, float y)
{
	return sqrt(x*x + y*y);
}

vector<float>
Gradient2(const vector<float>& D,
const Particle* p,
const int* dims)
{
	int x = p->x;
	int y = p->y;
	float dval1 = GetData2(D, x, y, dims[0], dims[1], (float)0);
	float dvalW = GetData2(D, x - 1, y, dims[0], dims[1], (float)dval1);
	float dvalE = GetData2(D, x + 1, y, dims[0], dims[1], (float)dval1);
	float dvalN = GetData2(D, x, y - 1, dims[0], dims[1], (float)dval1);
	float dvalS = GetData2(D, x, y + 1, dims[0], dims[1], (float)dval1);
	vector<float> g(2);
	g[0] = dvalE - dvalW;
	g[1] = dvalS - dvalN;
	float len = Length2(g[0], g[1]);
	g[0] /= len; g[1] /= len;
	return g;
}

vector<Particle*>
pickCore(vector<Particle*>& particles)
{
	vector<Particle*> core;
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
vector<vector<Particle*>>
clusterParticles(vector<Particle*>& particles, float thres)
{
	vector<Node<Particle*>*> nodes;
	for (int i = 0; i < particles.size(); ++i)
	{
		nodes.push_back(makeset(particles[i]));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		Particle* p = nodes[i]->key;
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			Particle* q = nodes[j]->key;
			if (Length2(p->x - q->x, p->y - q->y) <= thres)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<Particle*>*> rep = clusters(nodes);
	vector<vector<Particle*>> group(rep.size());
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
Cluster CParticle4D using disjoint set.
*/
vector<vector<CParticle4D>>
clusterParticles(vector<CParticle4D>& particles, float thres)
{
	vector<Node<CParticle4D>*> nodes;
	for (int i = 0; i < particles.size(); ++i)
	{
		nodes.push_back(makeset(particles[i]));
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CParticle4D p = nodes[i]->key;
		for (int j = i + 1; j < nodes.size(); ++j)
		{
			CParticle4D q = nodes[j]->key;
			if (Length2(p.m_X - q.m_X, p.m_Y - q.m_Y) <= thres)
			{
				merge(nodes[i], nodes[j]);
			}
		}
	}
	vector<Node<CParticle4D>*> rep = clusters(nodes);
	vector<vector<CParticle4D>> group(rep.size());
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
Check if there is a neighbor particle of the same generation that still has a room to grow. 
If so, then use the particle as a decendent of the core -> return true.
Else, this core is a legitimate one and should not be removed -> return false.
*/
bool
assignExit(vector<Particle*> core, //pass a copy
			vector<Particle*>& M, const int* dims)
{
	while (core.empty() == false)
	{
		vector<Particle*> core2;
		for (int i = 0; i < core.size(); ++i)
		{
			float mind = std::numeric_limits<float>::infinity();
			Particle* px = NULL;
			for (int n = 0; n < NumNeighborsHere; ++n)
			{
				CParticle4D p0(core[i]->x + XOffsetHere[n], core[i]->y + YOffsetHere[n]);
				Particle* q = GetVoxel(M, p0, (Particle*)NULL, dims);
				if (q && q->gen == core[i]->gen && q->descendents.empty() == false)
				{
					float d = Length2(core[i]->x - p0.m_X, core[i]->y - p0.m_Y);
					if (d < mind)
					{
						mind = d;
						px = q;
					}
				}
			}
			if (px == NULL)
			{
				core2.push_back(core[i]);
			}
			else
			{
				core[i]->descendents.insert(px);
				px->ascendents.insert(core[i]);
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
vector<Particle*>
removeIsolatedCore(vector<Particle*> core, //pass a copy
			int thres, vector<Particle*>& M, const int* dims)
{
	vector<vector<Particle*>> group = clusterParticles(core, 1.0f);
	//for each cluster, check if there is a neighbor that can be grown.
	//If so, assign its decentent to the neighbor (done all in assignExit).
	for (int g = 0; g < group.size(); ++g)
	{
		if (group[g].size() < thres)
		{
			if (assignExit(group[g], M, dims) == true) //there is an exit particle
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

vector<Particle*>
propagateParticles(vector<int>& S,
vector<unsigned char>& L,
vector<Particle*>& front,
int mincore,
const int* dims)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<Particle*> particles;
	vector<Particle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	vector<Particle*> Q;
	for (int i = 0; i < front.size(); ++i)
	{
		Q.push_back(front[i]);
		SetVoxel(mp, front[i], front[i], dims);
		particles.push_back(front[i]);
	}
	//int mincore = 5;
	int gen = 1;
	while (Q.empty() == false)
	{
		vector<Particle*> core;
		for (int i = 0; i < Q.size(); ++i)
		{
			SetVoxel(S, Q[i], gen, dims);
		}
		//G.push_back(Q); //separate per generation

		vector<Particle*> Q2;
		for (int i = 0; i<Q.size(); ++i)
		{
			for (int n = 0; n < NumNeighborsHere; ++n)
			{
				CParticle4D p0(Q[i]->x + XOffsetHere[n], Q[i]->y + YOffsetHere[n]);
				if (GetVoxel(L, p0, (unsigned char)0, dims) && GetVoxel(S, p0, 1, dims) == 0)
				{
					Particle* p = GetVoxel(mp, p0, (Particle*)NULL, dims);
					if (p == NULL)
					{
						p = factory.makeParticle(p0.m_X, p0.m_Y, gen);
						SetVoxel(mp, p0, p, dims);
						Q2.push_back(p);
						particles.push_back(p);
					}
					p->ascendents.insert(Q[i]);
					Q[i]->descendents.insert(p);
				}
			}
			if (Q[i]->descendents.empty())
			{
				core.push_back(Q[i]);
			}
		}
		core = removeIsolatedCore(core, mincore, mp, dims);

		Q = Q2;
		gen++;
	}
	return particles;
}

map<pair<int, int>, int> 
cutCores(map<pair<int, int>, int>& pairs, vector<int>& occupied, vector<int>& mask, const int* dims)
{
	map<pair<int, int>, int> cuts;
	vector<unsigned char> L(numberOfElements(2, dims));
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
	while (true)
	{
		vector<int> C(L.size(), 0);
		int nc = ConnectedComponentAnalysisBigger(C, L, NeighborhoodFour, (unsigned char)0, 2, dims);
		vector<Node<int>*> nodes;
		map<int, Node<int>*> nmap;
		for (set<int>::iterator it = labels.begin(); it != labels.end(); ++it)
		{
			Node<int>* n = makeset(*it);
			nodes.push_back(n);
			nmap[n->key] = n;
		}

		vector<set<int>> group(nc+1);
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
		for (map<pair<int, int>, int>::iterator it = pairs.begin(); it != pairs.end(); ++it)
		{
			pair<int, int> pr = it->first;
			if (cuts.find(pr) == cuts.end())
			{
				if (findset(nmap[pr.first]) != findset(nmap[pr.second]))
				{
					cuts[pr] = iter;
				}
			}
		}
		for (int i = 0; i < nodes.size(); ++i)
		{
			delete nodes[i];
		}
		
		vector<CParticle4D> rem;
		for (int i = 0; i < dims[1]; ++i)
		{
			for (int j = 0; j < dims[0]; ++j)
			{
				unsigned char val = GetData2(L, j, i, dims[0], dims[1], (unsigned char)0);
				int val2 = GetData2(mask, j, i, dims[0], dims[1], -1);
				if (val && val2 == -1)
				{
					if (GetData2(L, j - 1, i, dims[0], dims[1], (unsigned char)0) == 0 ||
						GetData2(L, j + 1, i, dims[0], dims[1], (unsigned char)0) == 0 ||
						GetData2(L, j, i - 1, dims[0], dims[1], (unsigned char)0) == 0 ||
						GetData2(L, j, i + 1, dims[0], dims[1], (unsigned char)0) == 0)
					{
						rem.push_back(CParticle4D(j, i));
					}

				}
			}
		}
		if (rem.empty()) break;
		for (int i = 0; i < rem.size(); ++i)
		{
			SetVoxel(L, rem[i], (unsigned char)0, dims);
		}
		iter++;
	}

	return cuts;
}


vector<vector<Particle*>>
groupCoreParticles(vector<Particle*>& particles, 
					vector<vector<Particle*>>& init,
					vector<int>& M, 
					float thres, const int* dims)
{
	vector<Particle*> pmap(numberOfElements(2, dims));
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->corelink = NULL;
		SetVoxel(pmap, particles[i], (Particle*)particles[i], dims);
		particles[i]->label = -1; //label for non-core
	}
	vector<int> occupied(numberOfElements(2, dims), -1);
	vector<Particle*> Q;
	for (int i = 0; i < init.size(); ++i)
	{
		for (int j = 0; j < init[i].size(); ++j)
		{
			SetVoxel(occupied, init[i][j], i, dims);
			init[i][j]->label = i;
			Q.push_back(init[i][j]);
		}
	}
	map<pair<int, int>, int> pairs;
	int m = 0;
	set<Particle*> done;
	while (Q.empty() == false)
	{
		m++;
		set<Particle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			Particle* p = Q[i];
			SetVoxel(M, p, m, dims);
			for (int n = 0; n < NumNeighborsFour; ++n)
			{
				CParticle4D qc(p->x + XOffsetFour[n], p->y + YOffsetFour[n]);
				Particle* q = GetVoxel(pmap, qc, (Particle*)NULL, dims);
				if (q != NULL)
				{
					int k = GetVoxel(occupied, qc, -1, dims);
					if (k < 0)
					{
						q->label = p->label;
						SetVoxel(occupied, qc, p->label, dims);
						Q2.insert(q);
					}
					else if (k != p->label)
					{
						pair<int, int> pr(Min(k, p->label), Max(k, p->label));
						if (pairs.find(pr) == pairs.end())
						{
							pairs[pr] = m;
						}
					}
				}
			}
		}
		Q.clear();
		Q.insert(Q.begin(), Q2.begin(), Q2.end());
	}

	vector<int> mask(numberOfElements(2, dims), -1);
	for (int i = 0; i < init.size(); ++i)
	{
		for (int j = 0; j < init[i].size(); ++j)
		{
			SetVoxel(mask, init[i][j], i, dims); //used later in erosion
			init[i][j]->label = i;
			Q.push_back(init[i][j]);
		}
	}
	map<pair<int, int>, int> cuts = cutCores(pairs, occupied, mask, dims);
	//now propagate back from each core to create a tree, and merge/cut the edge based on the cut-cost.
	vector<Node<int>*> nodes;
	for (int i = 0; i < init.size(); ++i)
	{
		nodes.push_back(makeset(i));
	}
	for (map<pair<int, int>, int>::iterator it = pairs.begin(); it != pairs.end(); ++it)
	{
		pair<int, int> pr = it->first;
		int dist = it->second;
		int cut = cuts[pr];
		if ((float)dist / cut < thres)
		{
			merge(nodes[pr.first], nodes[pr.second]);
		}
		printf("%d - %d: %d, %d, %f\n", pr.first, pr.second, cut, dist, (float)dist / cut);
	}
	vector<Node<int>*> reps = clusters(nodes);
	vector<vector<Particle*>> group(reps.size());
	{
		for (int i = 0; i < init.size(); ++i)
		{
			int k = distance(reps.begin(), find(reps.begin(), reps.end(), findset(nodes[i])));
			for (int j = 0; j < init[i].size(); ++j)
			{
				group[k].push_back(init[i][j]);
			}
		}
	}
	//reset label
	for (int i = 0; i < particles.size(); ++i)
	{
		particles[i]->label = 0;
	}
	//clean up 
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return group;
}

set<Particle*> trace(Particle* p)
{
	set<Particle*> S;
	set<Particle*> Q;
	Q.insert(p);
	while (Q.empty() == false)
	{
		set<Particle*> next;
		for (set<Particle*>::iterator it = Q.begin(); it != Q.end(); ++it)
		{
			Particle* q = *it;
			if (q->label == 0) q->label = -2;
			S.insert(q);
			for (set<Particle*>::iterator it = q->ascendents.begin(); it != q->ascendents.end(); ++it)
			{
				Particle* a = *it;
				next.insert(a);
				if (a->label == 0) a->label = q->label;
				//else if(a->label != q->label)  a->label = -1; //conflict.
			}
		}
		Q = next;
	}
	return S;
}

set<Particle*>
colorParticles(vector<vector<Particle*>>& cores, int offset)
{
	set<Particle*> S;
	for (int i = 0; i < cores.size(); ++i)
	{
		int label = i + 1;
		for (int j = 0; j < cores[i].size(); ++j)
		{
			Particle* p = cores[i][j];
			p->label = label + offset;
			set<Particle*> s = trace(p);
			S.insert(s.begin(), s.end());
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

	float thres = 2.0f;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(thres, prhs[1], classMode);
	}
	int minsize = 5;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minsize, prhs[2], classMode);
	}

	ParticleFactory& factory = ParticleFactory::getInstance();
	int nvoxels = numberOfElements(ndimL, dimsL);

	//vector<Particle*> front;
	//vector<Particle*> core;
	//set<Particle*> colored;
	map<Particle*, int> initgroup; //initial core group
	map<Particle*, int> finalgroup; //final core  group
	vector<int> S(nvoxels, 0);
	vector<int> S2(nvoxels, 0);

	vector<Particle*> front = ExtactSurface(L, dimsL);
	vector<Particle*> particles = propagateParticles(S, L, front, minsize, dimsL);
	vector<Particle*> core = pickCore(particles);
	vector<vector<Particle*>> init = clusterParticles(core, 1.1f); //initial clusters
	vector<vector<Particle*>> coregroup = groupCoreParticles(particles, init, S2, thres, dimsL);
	set<Particle*> colored = colorParticles(coregroup, 0);

	for (int i = 0; i < init.size(); ++i)
	{
		for (int j = 0; j < init[i].size(); ++j)
		{
			initgroup[init[i][j]] = i;
		}
	}
	for (int i = 0; i < coregroup.size(); ++i)
	{
		for (int j = 0; j < coregroup[i].size(); ++j)
		{
			finalgroup[coregroup[i][j]] = i;
		}
	}

	if (nlhs >= 1)
	{
		int dims[] = { core.size(), 7 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < core.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], core[i]->x);
			SetData2(F, i, 1, dims[0], dims[1], core[i]->y);
			SetData2(F, i, 2, dims[0], dims[1], core[i]->label);
			SetData2(F, i, 3, dims[0], dims[1], core[i]->gen);
			SetData2(F, i, 4, dims[0], dims[1], core[i]->id);
			SetData2(F, i, 5, dims[0], dims[1], initgroup[core[i]]);
			SetData2(F, i, 6, dims[0], dims[1], finalgroup[core[i]]);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { colored.size(), 5 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (set<Particle*>::iterator it = colored.begin(); it != colored.end(); ++it, ++i)
		{
			Particle* p = *it;
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->label);
			SetData2(F, i, 3, dims[0], dims[1], p->gen);
			SetData2(F, i, 4, dims[0], dims[1], p->id);
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
	mexUnlock();
	factory.clean();
}

