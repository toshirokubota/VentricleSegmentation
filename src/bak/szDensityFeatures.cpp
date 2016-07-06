
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>
using namespace std;

#include <szMiscOperations.h>
#include <szMexUtility.h>	
#include <szConnectedComponent.h>
#include <szCutOutNodule.h>
#include <szLocalExtrema.h>
#include <szMexUtilityTemplate.h>
#include <szEigen.h>
#include <szFeatureClass.h>

inline int
NumDensityFeatures()
{
	const int NumFeatures = 7;
	return NumFeatures;
}

vector<CFeature>
computeDensityFeatures(const vector<unsigned short>& A, //image
					   const vector<unsigned char>& B, //nodule segmentation
					   int ndim,
					   const int* dims) 
{
	int i;
	int nvoxels=numberOfElements(ndim,dims);
	vector<unsigned short> values;
	int vols=0; //# voxels in non-convex (nC) volume
	double sum = 0;
	double sum2 = 0;
	int cntBelow600 = 0;
	for(i=0; i<nvoxels; ++i) 
	{
		if(B[i]) {
			unsigned short val = A[i];
			sum += val;
			sum2 += val*val;
			values.push_back(val);
			if(val <= 600)
			{
				cntBelow600++;
			}
		}
	}
	vector<CFeature> features(NumDensityFeatures());
	if(values.size())
	{
		double mean = sum/values.size();
		double var = sum2/values.size() - mean*mean;
		sort(values.begin(), values.end());
		double median = values[values.size()/2];
		double q1 = values[(int)(values.size()*0.25)];
		double q3 = values[(int)(values.size()*0.75)];
		double p90 = values[(int)(values.size()*0.9)];
		features[0].Value = (double)cntBelow600/(double)values.size();
		features[0].Name = "SegCntBelow600";
		features[1].Value = mean/1000.0;
		features[1].Name = "SegMeanGray";
		features[2].Value = sqrt(var)/100.0;
		features[2].Name = "SegStdGray";
		features[3].Value = median/1000.0;
		features[3].Name = "SegMedianGray";
		features[4].Value = q1/1000.0;
		features[4].Name = "SegFirstQGray";
		features[5].Value = q3/1000.0;
		features[5].Name = "SegThirdQGray";
		features[6].Value = p90/1000.0;
		features[6].Name = "Seg90PGray";
	}
	return features;
}

