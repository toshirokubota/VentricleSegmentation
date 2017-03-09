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
	}
	int x;
	int y;
	int label;
	int gen; //generation
	int id;
	set<Particle*> ascendents;
	set<Particle*> descendents;
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
	vector<Particle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	for (int i = 0; i < front.size(); ++i)
	{
		SetVoxel(mp, front[i], front[i], dims);
	}
	map<Particle*, float> gradx;
	map<Particle*, float> grady;
	vector<Particle*> front2;
	/*
	find a direction for the growth. Create more particles if needed.
	*/
	for (int i = 0; i < front.size(); ++i)
	{
		vector<CParticle4D> region;
		for (int n = 0; n < NumNeighborsHere; ++n)
		{
			CParticle4D p(front[i]->x + XOffsetHere[n], front[i]->y + YOffsetHere[n]);
			if (GetVoxel(L, p, (unsigned char)0, dims) && GetVoxel(mp, p, (Particle*)NULL, dims) == NULL)
			{
				region.push_back(p);
			}
		}
		vector<vector<CParticle4D>> g = clusterParticles(region, 1.5f);
		if (g.size() == 0) //isolated point
		{
			gradx[front[i]] = 0;
			grady[front[i]] = 0;
		}
		else
		{
			vector<Particle*> ps(1, front[i]);
			for (int j = 1; j < g.size(); ++j)
			{
				//create additional particles for disjoint regions
				ps.push_back(factory.makeParticle(front[i]->x, front[i]->y, front[i]->gen, front[i]->label));
			}
			for (int j = 0; j < ps.size(); ++j)
			{
				vector<float> dir(2, 0.0f);
				for (int k = 0; k < g[j].size(); ++k)
				{
					dir[0] += (g[j][k].m_X - ps[j]->x);
					dir[1] += (g[j][k].m_Y - ps[j]->y);
				}
				dir[0] /= (float)g[j].size();
				dir[1] /= (float)g[j].size();
				gradx[ps[j]] = dir[0];
				grady[ps[j]] = dir[1];
				front2.push_back(ps[j]);
			}
		}
	}

	vector<Particle*> cores;
	vector<Particle*> Q;
	for (int i = 0; i < front.size(); ++i)
	{
		Q.push_back(front[i]);
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
		vector<Particle*> Q2;
		for (int i = 0; i<Q.size(); ++i)
		{
			float dx = gradx[Q[i]];
			float dy = grady[Q[i]];
			for (int n = 0; n < NumNeighborsHere; ++n)
			{
				CParticle4D p0(Q[i]->x + XOffsetHere[n], Q[i]->y + YOffsetHere[n]);
				if (GetVoxel(L, p0, (unsigned char)0, dims) && GetVoxel(S, p0, 1, dims) == 0)
				{
					if (XOffsetHere[n] * dx + YOffsetHere[n] * dy > 0)
					{
						Particle* p = GetVoxel(mp, p0, (Particle*)NULL, dims);
						if (p == NULL)
						{
							p = factory.makeParticle(p0.m_X, p0.m_Y, gen);
							SetVoxel(mp, p0, p, dims);
							Q2.push_back(p);
						}
						p->ascendents.insert(Q[i]);
						Q[i]->descendents.insert(p);
					}
				}
			}
			if (Q[i]->descendents.empty())
			{
				core.push_back(Q[i]);
			}
		}
		core = removeIsolatedCore(core, mincore, mp, dims);
		cores.insert(cores.end(), core.begin(), core.end());

		for (int i = 0; i < Q2.size(); ++i)
		{
			float dir[2] = { 0.0f, 0.0f };
			Particle* q = Q2[i];
			for (set<Particle*>::iterator it = q->ascendents.begin(); it != q->ascendents.end(); ++it)
			{
				dir[0] += gradx[*it];
				dir[1] += grady[*it];
			}
			gradx[q] = dir[0] / q->ascendents.size();
			grady[q] = dir[1] / q->ascendents.size();
		}

		Q = Q2;
		gen++;
	}
	return cores;
}

/*
vector<Particle*>
propagateParticles(vector<int>& S,
vector<unsigned char>& L,
vector<Particle*>& front,
int mincore,
const int* dims)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<Particle*> cores;
	vector<Particle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	vector<Particle*> Q;
	for (int i = 0; i < front.size(); ++i)
	{
		Q.push_back(front[i]);
		SetVoxel(mp, front[i], front[i], dims);
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
		cores.insert(cores.end(), core.begin(), core.end());

		Q = Q2;
		gen++;
	}
	return cores;
}
*/

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
colorParticles(vector<vector<Particle*>>& cores)
{
	set<Particle*> S;
	for (int i = 0; i < cores.size(); ++i)
	{
		int label = i + 1;
		for (int j = 0; j < cores[i].size(); ++j)
		{
			Particle* p = cores[i][j];
			p->label = label;
			set<Particle*> s = trace(p);
			S.insert(s.begin(), s.end());
		}
	}
	return S;
}

/*
This function groun core particles based on the ancestry trees.
For each core, it goes up its tree for NGEN generations.
Two cores that share a ancestor are grouped together.
*/
vector<vector<Particle*>>
groupCoreParticles(vector<Particle*>& core, int ngen)
{
	vector<set<Particle*>> vancestors;
	for (int i = 0; i < core.size(); ++i)
	{
		set<Particle*> aset;
		set<Particle*> current;
		current.insert(core[i]);
		for (int j = 0; j < ngen; ++j)
		{
			set<Particle*> next;
			for (set<Particle*>::iterator it = current.begin(); it != current.end(); ++it)
			{
				Particle* q = *it;
				for (set<Particle*>::iterator jt = q->ascendents.begin(); jt != q->ascendents.end(); ++jt)
				{
					if (aset.find(*jt) == aset.end())
					{
						next.insert(*jt);
						aset.insert(*jt);
					}
				}
			}
			current = next;
		}
		vancestors.push_back(aset);
	}
	vector<Node<Particle*>*> nodes;
	for (int i = 0; i < core.size(); ++i)
	{
		nodes.push_back(makeset(core[i]));
	}
	for (int i = 0; i < core.size(); ++i)
	{
		for (int j = i + 1; j < core.size(); ++j)
		{
			if (findset(nodes[i]) != findset(nodes[j]))
			{
				vector<Particle*> vint;
				set_intersection(vancestors[i].begin(), vancestors[i].end(), vancestors[j].begin(), vancestors[j].end(), back_inserter(vint));
				if (vint.empty() == false)
				{
					merge(nodes[i], nodes[j]);
				}
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

	vector<Particle*> front = ExtactSurface(L, dimsL);
	vector<int> S(nvoxels, 0);
	vector<Particle*> core = propagateParticles(S, L, front, minsize, dimsL);
	//vector<vector<Particle*>> coregroup = clusterParticles(core, thres);
	vector<vector<Particle*>> coregroup = groupCoreParticles(core, thres);
	set<Particle*> colored = colorParticles(coregroup);

	if (nlhs >= 1)
	{
		int dims[] = { core.size(), 5 };
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < core.size(); ++i)
		{
			SetData2(F, i, 0, dims[0], dims[1], core[i]->x);
			SetData2(F, i, 1, dims[0], dims[1], core[i]->y);
			SetData2(F, i, 2, dims[0], dims[1], core[i]->label);
			SetData2(F, i, 3, dims[0], dims[1], core[i]->gen);
			SetData2(F, i, 4, dims[0], dims[1], core[i]->id);
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
	mexUnlock();
	factory.clean();
}

