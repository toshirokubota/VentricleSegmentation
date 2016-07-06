#ifndef	___SZ_PBS_FEATURES_H___
#define ___SZ_PBS_FEATURES_H___

#include <vector>
using namespace std;

string GetPBSFeatureName1(const char* name0);

string GetPBSFeatureName2(const char* name0);

bool
PBS_IndexVectorWrapperShell(vector<int>& S, 
						   const vector<float>& D, 
						   int nshells);

bool
PBS_IndexVectorWrapperWedge(vector<int>& S, //index volume
							const vector<unsigned char>& B, //segmentation
							double x, double y, double z, //segmentation center
							int nwedges,
							int nbins,
							const int* dims);

bool
PBS_IndexVector(vector<int>& S, 
				const vector<float>& Dsq, 
				int nshells);

template<class Item>
vector<CFeature>
PBS_ShellFeatures(const vector<Item>& A,
				  const vector<int>& D,
				  int numBins,
				  const char* name,
				  const int* dims);


#include <../../src/NoduleSegmenationSupport/szPBSFeaturesTemplate.cpp>

#endif /* ___SZ_PBS_FEATURES_H___ */