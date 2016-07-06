#include <Kruskal.h>

vector<struct_edge>
kruskal(vector<struct_edge>& e, int numNodes)
{
	sort(e.begin(), e.end());
	CDisjointSet set(numNodes);
	int i;
	for(i=0; i<numNodes; ++i)
	{
		set.makeset(i);
	}

	vector<struct_edge> u;
	for(int i=0; i<e.size() && u.size() < numNodes - 1; ++i)
	{
		int s1 = set.findset(e[i].v1);
		int s2 = set.findset(e[i].v2);;
		if(s1 >= 0 && s2>=0 && s1 != s2)
		{
			set.merge(e[i].v1, e[i].v2);
			u.push_back(e[i]);
		}
	}
	return u;
}

