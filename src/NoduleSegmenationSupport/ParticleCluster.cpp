#include <mex.h>
#include <ParticleCluster.h>
#include <szMiscOperations.h>
#include <cmath>
#include <algorithm>
using namespace std;
#include <DisjointSet.h>
#include <Kruskal.h>
#include <szMexUtility.h>

float 
DistanceBetweenParticles(const CParticle& p1, const CParticle& p2)
{
	float d = (p1.m_X-p2.m_X)*(p1.m_X-p2.m_X) + (p1.m_Y-p2.m_Y)*(p1.m_Y-p2.m_Y) + (p1.m_Z-p2.m_Z)*(p1.m_Z-p2.m_Z);
	return sqrt(d);
}

vector<struct_edge>
BreakTree(const vector<struct_edge>& MST)
{
	vector<float> vweights(MST.size());
	for(int i=0; i<MST.size(); ++i)
	{
		vweights[i] = MST[i].weight;
	}

	sort(vweights.begin(), vweights.end());
	float q1 = vweights[vweights.size()/4];
	float q3 = vweights[(3*vweights.size())/4];
	float cutoff = q3 + 1.5*(q3-q1);
	//cutoff = Max(cutoff, 3.0);
	printf("BreakTree: q1=%f, q3=%f, cutoff=%f.\n", q1, q3, cutoff);

	vector<struct_edge> edges;
	int cnt = 0;
	for(int i=0; i<MST.size(); ++i)
	{
		if(MST[i].weight <= cutoff)
		{
			edges.push_back(MST[i]);
			cnt++;
		}
	}
	printf("BreakTree: out of %d edges, %d of them are kept.\n", MST.size(), cnt);
	return edges;
}


vector<int>
ClusterParticles(const vector<CParticle>& particles)
{
	int N = particles.size();
	vector<struct_edge> vedges;
	for(int i=0; i<N; ++i)
	{
		for(int j=i+1; j<N; ++j)
		{
			vedges.push_back(struct_edge(i, j, DistanceBetweenParticles(particles[i], particles[j])));
		}
	}

	vector<int> vlabels(N);
	vector<struct_edge> MST = kruskal(vedges, N);

	vector<struct_edge> vedges2 = BreakTree(MST);

	CDisjointSet set(particles.size());
	for(int i=0; i<N; ++i)
	{
		set.makeset(i);
	}
	for(int i=0; i<vedges2.size(); ++i)
	{
		set.merge(vedges2[i].v1, vedges2[i].v2);
	}

	for(int i=0; i<N; ++i)
	{
		vlabels[i] = set.findset(i);
	}

	return vlabels;
}


