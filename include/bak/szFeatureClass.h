#ifndef ___FEATURE_CLASS_H___
#define ___FEATURE_CLASS_H___
#ifdef MEX_DLL
#include <mex.h>
#endif
#include <ostream>
#include <vector>
using namespace std;

class CFeature {
public:
	CFeature(const string& name = "Unkown", double val = 0)
	{
		Value = val;
		Name = name;
	}
	CFeature(const CFeature& feature)
	{
		Value = feature.Value;
		Name = feature.Name;
	}
	CFeature& operator =(const CFeature& feature)
	{
		Value = feature.Value;
		Name = feature.Name;
		return *this;
	}
	double Value;
	string Name;
};

ostream& operator << (ostream& out, const CFeature& f);

#ifdef MEX_DLL
mxArray*
StoreFeatures(const vector<CFeature>& features);
#endif

#endif /* ___FEATURE_CLASS_H___ */