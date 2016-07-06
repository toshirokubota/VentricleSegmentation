#ifndef ___SZ_CORRELATION_FEATURES_H___
#define ___SZ_CORRELATION_FEATURES_H___

#include <vector>
#include <string>
using namespace std;
#include <szFeatureClass.h>

inline string
GetCorrelationFeatureName1(const char* name0)
{
	string name = "VectorAngle";
	name += name0;
	return name;
	string name2 = "Correlation";
	name2 += name;
}

inline string
GetCorrelationFeatureName2(const char* name0)
{
	string name = "Correlation";
	name += name0;
	return name;
}

template<class Item1, class Item2>
vector<CFeature>
computeCorrelationFeatures(const vector<Item1>& A, //image 1
						   const vector<Item2>& B,  //image 2
						   const vector<unsigned char>& M, //mask 
						   const char* name,
						   const int* dims);

#include <../../src/NoduleSegmenationSupport/szCorrelationFeaturesTemplate.cpp>

#endif /* ___SZ_CORRELATION_FEATURES_H___ */