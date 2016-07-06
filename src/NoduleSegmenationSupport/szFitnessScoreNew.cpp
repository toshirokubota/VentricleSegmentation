
//#include <mex.h>
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>
using namespace std;

#include <szCutOutNodule.h>
#include <szSizeFeaturesV3.h>
#include <szDensityFeatures.h>
#include <szGeometryFeatures.h>
#include <szTextureFeatures.h>

int
NumberOfFeatures()
{
	return NumSizeFeatures() + NumDensityFeatures() + NumGeometryFeatures() + NumTextureFeatures();
}

vector<double>
fitnessScoreNew(const vector<unsigned short>& A, //image
				const vector<unsigned char>& B, //nodule segmentation
				const vector<unsigned char>& L, //foreground segmentation
				const vector<float>& D,         //distance map
				const vector<float>& Sp,         //sphericity
				const vector<unsigned char>& P,  //sphericity local maxima
				const vector<int>& vseeds,      //seed positions
				const vector<float>& v,         //voxel size
				int ndim,
				const int* dims) 
{
	vector<double> vFeatures;

	//Size features
	//printf("Computing size features...\r");
	vector<double> size = computeSizeFeaturesV3(B, v, ndim, dims);
	vFeatures.insert(vFeatures.end(), size.begin(), size.end());
	//printf("Computing size features... Done.\n");

	//Density features
	//printf("Computing density features...\r");
	vector<double> density = computeDensityFeatures(A, B, ndim, dims);
	vFeatures.insert(vFeatures.end(), density.begin(), density.end());
	//printf("Computing density features... Done.\n");

	//Geometry features
	//printf("Computing geometry features...\r");
	//vector<unsigned char> L2 = L;
	//selectCluster(L2, vseeds[0], vseeds[1], vseeds[2], dims);
	vector<double> geom = computeGeometryFeatures(B, L, D, vseeds, v, ndim, dims);
	vFeatures.insert(vFeatures.end(), geom.begin(), geom.end());
	//printf("Computing geometry features... Done.\n");

	//Texture features
	//printf("Computing texture features...\r");
	vector<double> texture = computeTextureFeatures(A, B, ndim, dims);
	vFeatures.insert(vFeatures.end(), texture.begin(), texture.end());
	//printf("Computing texture features... Done.\n");

	return vFeatures;
}
