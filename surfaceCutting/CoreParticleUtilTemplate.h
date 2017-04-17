#ifndef __CORE_PARTICLE_UTIL_TEMPLATE_H__
#define __CORE_PARTICLE_UTIL_TEMPLATE_H__

#include <CoreParticle.h>
#include <mex.h>

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

//collect all possible K items from a vector of N items.
template<class T>
vector<vector<T>>
combinations(vector<T>& items, int k)
{
	if (k <= 0 || items.empty()) return vector<vector<T>>(); //empty
	else if (k == items.size())
	{
		return vector<vector<T>>(1, items);
	}
	else if (k == 1)
	{
		vector<vector<T>> result;
		for (int i = 0; i < items.size(); ++i)
		{
			result.push_back(vector<T>(1, items[i]));
		}
		return result;
	}
	else
	{
		vector<T> items2;
		items2.insert(items2.end(), items.begin(), items.end() - 1);
		vector<vector<T>> result = combinations(items2, k);
		vector<vector<T>> result2 = combinations(items2, k - 1);
		for (int i = 0; i < result2.size(); ++i)
		{
			result2[i].push_back(items[items.size() - 1]);
		}
		result.insert(result.end(), result2.begin(), result2.end());
		return result;
	}
}

#endif /* __CORE_PARTICLE_UTIL_TEMPLATE_H__ */