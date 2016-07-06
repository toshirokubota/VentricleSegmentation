#ifndef	___SZ_PBS_FEATURES_H___
#define ___SZ_PBS_FEATURES_H___

#include <vector>
using namespace std;

inline string
GetPBSFeatureName1(const char* name0)
{
	string name = "PBSDif";
	name += name0;
	return name;
}

inline string
GetPBSFeatureName2(const char* name0)
{
	string name = "PBS";
	name += name0;
	return name;
}

inline int
NumPBSFeatures();

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

/*vector<CFeature>
computePBSFeatures(const vector<unsigned short>& A, //gray scale data
				   const vector<unsigned char>& B,	//nodule segmentation
				   const vector<float>& Ds,			//distance map within the segmentation
				   int x, int y, int z,				//seed location
				   const int* dims);*/

#endif /* ___SZ_PBS_FEATURES_H___ */