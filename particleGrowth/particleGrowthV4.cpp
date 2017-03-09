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
	CoreParticle(int x0 = 0, int y0 = 0, int gen0 = 0, int lb = 0, int id0 = 0)
	{
		x = x0;
		y = y0;
		gen = gen0;
		label = lb;
		id = id0;
		//corelink = NULL;
	}
	int x;
	int y;
	int label;
	int gen; //generation
	int id;
	set<CoreParticle*> ascendents;
	set<CoreParticle*> descendents;
	//CoreParticle* corelink; //used to link cores

};

struct ParticleFactory
{
public:
	static ParticleFactory& getInstance()
	{
		static ParticleFactory instance;
		return instance;
	}
	CoreParticle* makeParticle(int x, int y, int gen, int lb = 0)
	{
		CoreParticle* particle = new CoreParticle(x, y, gen, lb, _id++);
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

struct cluster_pair;
struct core_cluster
{
	core_cluster(int id0)
	{
		id = id0;
	}
	vector<CoreParticle*> particles;
	vector<cluster_pair*> pairs;
	int size() { return particles.size(); }
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
		return (float)length / (float)cut * Min(u->size(), v->size()); // / Max(u->size(), v->size());
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

core_cluster*
merge_cores(core_cluster* a, core_cluster* b)
{
	CoreClusterFactory& factory = CoreClusterFactory::getInstance();
	core_cluster* c = factory.makeCoreCluster();
	c->particles = a->particles;
	c->particles.insert(c->particles.end(), b->particles.begin(), b->particles.end());
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
const CoreParticle* p,
T value,
const int* dims)
{
	return SetData2(A, p->x, p->y, dims[0], dims[1], value);
}

template<class T>
T
GetVoxel(const vector<T>& A,
const CoreParticle* p,
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

vector<CoreParticle*>
ExtactSurface(const vector<unsigned char>& L, 
			  const int* dims)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<CoreParticle*> front;
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
			if (Length2(p->x - q->x, p->y - q->y) <= thres)
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
			vector<CoreParticle*>& M, const int* dims)
{
	while (core.empty() == false)
	{
		vector<CoreParticle*> core2;
		for (int i = 0; i < core.size(); ++i)
		{
			float mind = std::numeric_limits<float>::infinity();
			CoreParticle* px = NULL;
			for (int n = 0; n < NumNeighborsHere; ++n)
			{
				CParticle4D p0(core[i]->x + XOffsetHere[n], core[i]->y + YOffsetHere[n]);
				CoreParticle* q = GetVoxel(M, p0, (CoreParticle*)NULL, dims);
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
vector<CoreParticle*>
removeIsolatedCore(vector<CoreParticle*> core, //pass a copy
int thres, vector<CoreParticle*>& M, const int* dims)
{
	vector<vector<CoreParticle*>> group = clusterParticles(core, 1.0f);
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

/*
from each surface particle, move towrad the center and establish ascendent/descendent relation.
*/
vector<CoreParticle*>
propagateParticles(vector<int>& S,
vector<unsigned char>& L,
vector<CoreParticle*>& front,
int mincore,
const int* dims)
{
	ParticleFactory& factory = ParticleFactory::getInstance();
	vector<CoreParticle*> particles;
	vector<CoreParticle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	vector<CoreParticle*> Q;
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
		vector<CoreParticle*> core;
		for (int i = 0; i < Q.size(); ++i)
		{
			SetVoxel(S, Q[i], gen, dims);
		}
		//G.push_back(Q); //separate per generation

		vector<CoreParticle*> Q2;
		for (int i = 0; i<Q.size(); ++i)
		{
			for (int n = 0; n < NumNeighborsHere; ++n)
			{
				CParticle4D p0(Q[i]->x + XOffsetHere[n], Q[i]->y + YOffsetHere[n]);
				if (GetVoxel(L, p0, (unsigned char)0, dims) && GetVoxel(S, p0, 1, dims) == 0)
				{
					CoreParticle* p = GetVoxel(mp, p0, (CoreParticle*)NULL, dims);
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

/*
For each adjacent pair of clusters, find the cut cost, i.e. the number of erosion operations to separate the two.
*/
map<cluster_pair*, int> 
clusterCutCost(vector<cluster_pair*>& pairs, 
				vector<core_cluster*>& init, //initial groups of core particles
				vector<int>& occupied, 
				const int* dims)
{
	vector<int> mask(numberOfElements(2, dims), -1);
	for (int i = 0; i < init.size(); ++i)
	{
		for (int j = 0; j < init[i]->size(); ++j)
		{
			SetVoxel(mask, init[i]->particles[j], i, dims); //used later in erosion
			init[i]->particles[j]->label = i;
		}
	}

	map<cluster_pair*, int> cuts;
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

		//for each cluster, check for the source (using occupied)
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
		
		//perform erosion - but do not erode initial cores.
		vector<CParticle4D> rem;
		for (int i = 0; i < dims[1]; ++i)
		{
			for (int j = 0; j < dims[0]; ++j)
			{
				unsigned char val = GetData2(L, j, i, dims[0], dims[1], (unsigned char)0);
				int val2 = GetData2(mask, j, i, dims[0], dims[1], -1); //check if it is an initial core
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

map<cluster_pair*, int>
clusterSeparationLength(const vector<CoreParticle*>& particles, //all particles
						const vector<core_cluster*>& init, //initial groups of core particles
						vector<int>& occupied, //mark the region growth and its origin
						const int* dims)
{
	vector<int> M(numberOfElements(2, dims), 0); //to store region-growth progress
	vector<CoreParticle*> pmap(numberOfElements(2, dims));
	for (int i = 0; i < particles.size(); ++i)
	{
		//particles[i]->corelink = NULL;
		SetVoxel(pmap, particles[i], (CoreParticle*)particles[i], dims);
		particles[i]->label = -1; //label for non-core
	}
	vector<CoreParticle*> Q;
	map<int, core_cluster*> cmap;
	for (int i = 0; i < init.size(); ++i)
	{
		cmap[i] = init[i];
		for (int j = 0; j < init[i]->size(); ++j)
		{
			SetVoxel(occupied, init[i]->particles[j], i, dims);
			init[i]->particles[j]->label = init[i]->id;
			Q.push_back(init[i]->particles[j]);
		}
	}
	CoreClusterFactory& factory = CoreClusterFactory::getInstance();
	vector<cluster_pair*> pairs;
	map<cluster_pair*, int> separation; //maps adjacent pairs and its separation (based on the time it took for them to meet).
	int m = 0;
	while (Q.empty() == false)
	{
		m++;
		set<CoreParticle*> Q2;
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			SetVoxel(M, p, m, dims);
			for (int n = 0; n < NumNeighborsFour; ++n)
			{
				CParticle4D qc(p->x + XOffsetFour[n], p->y + YOffsetFour[n]);
				CoreParticle* q = GetVoxel(pmap, qc, (CoreParticle*)NULL, dims);
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
		Q.clear();
		Q.insert(Q.begin(), Q2.begin(), Q2.end());
	}
	return separation;
}
/*
This routine first establish a graph by region-growing from each core - when two regions originated from different cores touch,
then the cores are considered adjacent.
It then perform grouping using the graph -- if the edge strengh is less than thres, then the nodes are put into the same group.
*/
vector<core_cluster*>
groupCoreParticles(vector<CoreParticle*>& particles, //all particles
					vector<core_cluster*>& init, //initial groups of core particles
					int nclusters, //desired number of clusters 
					const int* dims)
{
	vector<int> occupied(numberOfElements(2, dims), -1);
	map<cluster_pair*, int> separation = clusterSeparationLength(particles, init, occupied, dims);
	vector<cluster_pair*> core_pairs; //pairs of adjacent core clusters
	for (map<cluster_pair*, int>::iterator it = separation.begin(); it != separation.end(); ++it)
	{
		core_pairs.push_back(it->first);
	}
	map<cluster_pair*, int> cuts = clusterCutCost(core_pairs, init, occupied, dims);

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
	}

	while (clusters.size() > nclusters)
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
		if (mincost < std::numeric_limits<float>::infinity())
		{
			cluster_pair* pr = clusters[choice.first]->pairs[choice.second];
			printf("Merging: %d (%d) and %d (%d) - %d, %d, %f\n",
				pr->u->id, pr->u->size(), pr->v->id, pr->v->size(), (int)pr->length, (int)pr->cut, pr->weight());
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

set<CoreParticle*> trace(CoreParticle* p)
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
				if (a->label == 0) a->label = q->label;
				//else if(a->label != q->label)  a->label = -1; //conflict.
			}
		}
		Q = next;
	}
	return S;
}

set<CoreParticle*>
colorParticles(vector<core_cluster*>& cores, int offset)
{
	set<CoreParticle*> S;
	for (int i = 0; i < cores.size(); ++i)
	{
		int label = i + 1;
		for (int j = 0; j < cores[i]->size(); ++j)
		{
			CoreParticle* p = cores[i]->particles[j];
			p->label = label + offset;
			set<CoreParticle*> s = trace(p);
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

	int nclusters = 2;
	if (nrhs >= 2)
	{
		mxClassID classMode;
		ReadScalar(nclusters, prhs[1], classMode);
	}
	int minsize = 5;
	if (nrhs >= 3)
	{
		mxClassID classMode;
		ReadScalar(minsize, prhs[2], classMode);
	}

	ParticleFactory& factory = ParticleFactory::getInstance();
	int nvoxels = numberOfElements(ndimL, dimsL);

	vector<int> S(nvoxels, 0);

	vector<CoreParticle*> front = ExtactSurface(L, dimsL);
	vector<CoreParticle*> particles = propagateParticles(S, L, front, minsize, dimsL);
	vector<CoreParticle*> core = pickCore(particles);
	vector<core_cluster*> init = clusterCoreParticles(core, 1.1f); //initial clusters
	vector<core_cluster*> coregroup = groupCoreParticles(particles, init, nclusters, dimsL);
	set<CoreParticle*> colored = colorParticles(coregroup, 0);

	if (nlhs >= 1)
	{
		int count = 0;
		for (int i = 0; i < init.size(); ++i)
		{
			count += init[i]->size();
		}
		int dims[] = { count, 7 };
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
				SetData2(F, k, 2, dims[0], dims[1], p->label);
				SetData2(F, k, 3, dims[0], dims[1], p->gen);
				SetData2(F, k, 4, dims[0], dims[1], p->id);
				SetData2(F, k, 5, dims[0], dims[1], i);
				SetData2(F, k, 6, dims[0], dims[1], cmap[p]);
			}
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { colored.size(), 5 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (set<CoreParticle*>::iterator it = colored.begin(); it != colored.end(); ++it, ++i)
		{
			CoreParticle* p = *it;
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
		int dims[] = { core.size(), 5 };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (vector<CoreParticle*>::iterator it = core.begin(); it != core.end(); ++it, ++i)
		{
			CoreParticle* p = *it;
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->label);
			SetData2(F, i, 3, dims[0], dims[1], p->gen);
			SetData2(F, i, 4, dims[0], dims[1], p->id);
		}
		plhs[3] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	mexUnlock();
	factory.clean();
	CoreClusterFactory::getInstance().clean();
}

