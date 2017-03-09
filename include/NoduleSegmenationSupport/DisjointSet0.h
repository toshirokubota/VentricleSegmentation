#ifndef _DISJOINT_SET_H_
#define _DISJOINT_SET_H_

#include <vector>
#include <list>
using namespace std;

struct CDisjointSet
{
	CDisjointSet(int n=0)
	{
		vParent = vector<int>(n, -1);
		vRank = vector<int>(n, 0);
	};
	bool makeset(int i);
	bool merge(int i, int j);
	int findset(int i);

	vector<int> vParent;
	vector<int> vRank;
};

#endif /* _DISJOINT_SET_H_ */