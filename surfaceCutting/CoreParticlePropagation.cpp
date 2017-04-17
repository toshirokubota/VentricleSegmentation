#include <CoreParticlePropagation.h>
#include <CoreParticleUtil.h>
#include <CoreParticleUtilTemplate.h>
#include <Kruskal.h>
#include <Graph.h>
#include <GraphFactory.h>
#include <map>
using namespace std;

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

inline int COMPARE(int a, int b)
{
	//return a > b ? 1 : (a < b ? -1 : 0);
	return a == b ? 2 : (a < b ? 1 : 0);
}

int neighborScore(CoreParticle* p, CoreParticle* q)
{
	int x = COMPARE(p->x, q->x) + 1;
	int y = COMPARE(p->y, q->y) + 1;
	int z = COMPARE(p->z, q->z) + 1;
	int t = COMPARE(p->t, q->t) + 1;

	return x + 3 * y + 9 * z + 27 * t;
}

void
simplifyDescendents(vector<CoreParticle*>& particles)
{
	map<CoreParticle*, int> score;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (particles[i]->descendents.size() > 1)
		{
			score[particles[i]] = 0;
		}
		else
		{
			score[particles[i]] = 1;
		}
	}
	while (true)
	{
		bool bChanged = false;
		for (int i = 0; i < particles.size(); ++i)
		{
			CoreParticle* p = particles[i];
			if (p->descendents.size() > 1)
			{
				vector<pair<int, CoreParticle*>> toremove;
				for (set<CoreParticle*>::iterator it = p->descendents.begin(); it != p->descendents.end(); ++it)
				{
					CoreParticle* q = *it;
					for (set<CoreParticle*>::iterator jt = q->ascendents.begin(); jt != q->ascendents.end(); ++jt)
					{
						CoreParticle* r = *jt;
						if (r->descendents.size() == 1) // || neighborScore(p, q) < neighborScore(r, q))
						{
							toremove.push_back(pair<int, CoreParticle*>(0, q));
							break;
						}
						else if (neighborScore(p, q) < neighborScore(r, q))
						{
							toremove.push_back(pair<int, CoreParticle*>(neighborScore(p, q), q));
							break;
						}
					}
				}
				sort(toremove.begin(), toremove.end());
				for (int k = 0; k < Min(toremove.size(), p->descendents.size() - 1); ++k)
				{
					p->descendents.erase(toremove[k].second);
					(toremove[k].second)->ascendents.erase(p);
					bChanged = true;
				}
			}
		}
		if (bChanged == false) break;
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
		//p->value = 1; //TK!!!
		p->value = 0;
		if (surfaceParticle(p, ndim))
		{
			Q.push_back(p);
			p->value = 1;
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
					p->descendents.insert(q);
					q->ascendents.insert(p);
				}
			}
		}
		propagateDescendency(Q);

		simplifyDescendents(Q);
		for (int i = 0; i < Q.size(); ++i)
		{
			CoreParticle* p = Q[i];
			//distribute the value to its decendents
			if (p->descendents.empty() == false)
			{
				float amnt = p->value / p->descendents.size();
				for (set<CoreParticle*>::iterator it = p->descendents.begin(); it != p->descendents.end(); ++it)
				{
					CoreParticle* q = *it;
					q->value += amnt;
					//q->value = Min(q->value, amnt); //TK!!!
				}
				//p->value = 0.0f;
			}
		}

		Q.clear();
		Q.insert(Q.end(), Q2.begin(), Q2.end());

		gen++;
	}
	return particles;
}


/*
Connect medial axis particles as minimum spanning tree.
*/
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
			if (medialParticle(mp[i], ndim))
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
				{					CoreParticle* q = *it;
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
	
	return mst;
}
