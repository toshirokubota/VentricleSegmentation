
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>
using namespace std;

#include <szMiscOperations.h>
#include <szMexUtility.h>	
#include <szMexUtilityTemplate.h>
#include <szDistributionsTemplate.h>
#include <szHybridFeatures.h>
#include <szFeatureClass.h>

vector<CFeature>
hybridFeatures()
{
	int numFeatures = 1;
	vector<CFeature> features(numFeatures);
	features[0].Name = "PseudoDivergence";
	return features;
}

vector<CFeature>
computeHybridFeatures(const vector<unsigned short>& A, //image
					  const vector<unsigned char>& B, //nodule segmentation
					  const vector<float>& D, 
					  const int* dims)
{
	vector<CFeature> features = hybridFeatures();
	int i, j, k;
	//get the centroid of the segmentation
	vector<double> vcent = computeCentroid(B, 3, dims);
	float dval = GetData3(D, (int)vcent[0], (int)vcent[1], (int)vcent[2], dims[0], dims[1], dims[2], (float)0);

	//compute a divergence feature
	double dvg = 0;
	for(i=0; i<dims[2]; ++i) 
	{
		for(j=0; j<dims[1]; ++j) 
		{
			for(k=0; k<dims[0]; ++k) 
			{
				if(GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					unsigned short c = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (unsigned short)0);
					unsigned short w = GetData3(A, k-1, j, i, dims[0], dims[1], dims[2], c);
					unsigned short e = GetData3(A, k+1, j, i, dims[0], dims[1], dims[2], c);
					unsigned short n = GetData3(A, k, j-1, i, dims[0], dims[1], dims[2], c);
					unsigned short s = GetData3(A, k, j+1, i, dims[0], dims[1], dims[2], c);
					unsigned short f = GetData3(A, k, j, i-1, dims[0], dims[1], dims[2], c);
					unsigned short b = GetData3(A, k, j, i+1, dims[0], dims[1], dims[2], c);
					float gx = e - w;
					float gy = s - n;
					float gz = b - f;
					float dx = k - vcent[0];
					float dy = j - vcent[1];
					float dz = i - vcent[2];
					float dm = dx*dx + dy*dy + dz*dz;
					float wgt = exp(-dm/(2.*dval*dval));
					dvg -= wgt*(gx*dx + gy*dy + gz*dz);
				}
			}
		}
	}
	double beta = 0.00001;
	features[0].Value = atan(dvg*beta);

	return features;
}
