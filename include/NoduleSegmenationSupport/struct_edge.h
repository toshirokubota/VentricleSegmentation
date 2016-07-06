#ifndef _STRUCT_EDGE_H_
#define _STRUCT_EDGE_H_

#include <fstream>
using namespace std;

struct struct_edge
{
public:
	struct_edge(int v=0, int u=0, int w=0)
	{
		v1 = v;
		v2 = u;
		weight = w;
	}

	bool operator < (const struct_edge& e) const
	{
		if(weight < e.weight)
		{
			return true;
		}
		else
		{
			return false;
		}	
	}

	int v1;
	int v2;
	int weight;
};

#endif /* _STRUCT_EDGE_H_ */