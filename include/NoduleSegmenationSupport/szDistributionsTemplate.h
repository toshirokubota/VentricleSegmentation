#ifndef ___SZ_DISTRIBUTIONS_TEMPLATE_H___
#define ___SZ_DISTRIBUTIONS_TEMPLATE_H___

#include <vector>
using namespace std;
#include <szMexUtility.h>
#include <szFeatureClass.h>
#include <szMexUtilityTemplate.h>

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
double
JSD_Features(const vector<Item>& A,
			 const vector<Item>& B,
			 int numBins);

template<class Item1, class Item2>
float 
computeCorrelation(const vector<Item1>& va,
				   const vector<Item2>& vb,
				   int length = 0);

#include <../../src/NoduleSegmenationSupport/szDistributionsTemplate.cpp>

#endif /* ___SZ_DISTRIBUTIONS_TEMPLATE_H___ */