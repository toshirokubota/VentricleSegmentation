#include <DisjointSet.h>


bool
CDisjointSet::makeset(int i)
{
	if(i >= vParent.size())
	{
		vParent.resize(i + 1, -1);
		vRank.resize(i + 1, 0);
	}
	vParent[i] = i;
	return true;
}

int
CDisjointSet::findset(int n)
{
	if(n >= vParent.size())
	{
		return -1;
	}
	else if(n != vParent[n])
	{
		vParent[n] = findset(vParent[n]);
	}
	return vParent[n];
}

bool
CDisjointSet::merge(int n, int m)
{
	n = findset(n);
	m = findset(m);
	if(n>=0 && m>=0)
	{
		if(n == m)
		{
			return false;
		}
		else
		{
			if(vRank[m] > vRank[n])
			{
				vParent[n] = m;
			}
			else
			{
				vParent[m] = n;
				if(vRank[m] == vRank[n])
				{
					vRank[m]++;
				}
			}
			return true;
		}
	}
	else
	{
		return false;
	}
}
