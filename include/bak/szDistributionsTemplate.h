#ifndef ___SZ_DISTRIBUTIONS_TEMPLATE_H___
#define ___SZ_DISTRIBUTIONS_TEMPLATE_H___

#include <vector>
using namespace std;
#include <szMexUtility.h>
#include <szFeatureClass.h>

template<class Item>
vector<int>
computeFrequency(const vector<Item>& values,
				 const Item& minVal,
				 const Item& maxVal,
				 int numBins);

template<class Item>
vector<int>
computeFrequency(const vector<Item>& values,
				 const vector<Item>& vbins);

template<class Item>
vector<float>
computePDF(const vector<Item>& values,
		   const Item& minVal,
		   const Item& maxVal,
		   int numBins);

template<class Item>
vector<float>
computePDF(const vector<Item>& values,
		   const vector<Item>& vbins);

template<class Item>
vector<float>
computeCDF(const vector<Item>& values,
		   const Item& minVal,
		   const Item& maxVal,
		   int numBins);

template<class Item>
vector<float>
computeCDF(const vector<Item>& values,
		   const vector<Item>& vbins);

template<class Item>
bool
IndexVector(vector<int>& S, 
			const vector<Item>& B, 
			const vector<Item>& vedges);

template<class Item>
void
gradientMagnitude(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims);

template<class Item>
void
laplacianResponse(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims);

template<class Item>
void
gradientAngle(vector<float>& G,
			  const vector<Item>& A,
			  const vector<unsigned char>& S,
			  int x, int y, int z,
			  const int* dims);

template<class Item>
double
JSD_Features(const vector<Item>& A,
			 const vector<Item>& B,
			 int numBins);

template<class Item>
vector<CFeature>
PBS_ShellFeatures(const vector<Item>& A,
				  const vector<int>& D,
				  int numBins,
				  const char* name,
				  const int* dims);

template<class Item1, class Item2>
float 
computeCorrelation(const vector<Item1>& va,
				   const vector<Item2>& vb,
				   int length = 0);

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

#include <../src/szDistributionsTemplate.cpp>

#endif /* ___SZ_DISTRIBUTIONS_TEMPLATE_H___ */