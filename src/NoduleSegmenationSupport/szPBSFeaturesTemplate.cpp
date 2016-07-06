#include <vector>
using namespace std;
#include <szFeatureClass.h>
#include <szPBSFeatures.h>

inline
vector<CFeature>
pbsFeatures(const char* name)
{	
	int numFeatures = 2;
	vector<CFeature> features(numFeatures);
	features[0].Name = GetPBSFeatureName1(name);
	features[1].Name = GetPBSFeatureName2(name);

	return features;
}

template<class Item>
vector<CFeature>
PBS_ShellFeatures(const vector<Item>& A,
				  const vector<int>& D,
				  int numBins,
				  const char* name,
				  const int* dims)
{
	vector<CFeature> features = pbsFeatures(name);
	int i, j, k;
	unsigned short numskins = 0;
	for(i=0; i<D.size(); ++i)
	{
		if(D[i] > numskins)
		{
			numskins = D[i];
		}
	}

	vector<vector<Item>> vbuckets(numskins);
	vector<Item> allvals;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				int bin = GetData3(D, k, j, i, dims[0], dims[1], dims[2], (int)0);
				if(bin)
				{
					Item val = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (Item)0);
					vbuckets[bin - 1].push_back(val);
					allvals.push_back(val);
				}
			}
		}
	}
	if(!allvals.empty())
	{
		Item minVal = allvals[0];
		Item maxVal = allvals[0];
		for(i=1; i<allvals.size(); ++i)
		{
			minVal = Min(minVal, allvals[i]);
			maxVal = Max(maxVal, allvals[i]);
		}

		vector<float> pdfAll = computePDF(allvals, minVal, maxVal, numBins);

		vector<float> pdfUniform(pdfAll.size(), 1.0/pdfAll.size());
		float jsdAll = computeJSD(pdfAll, pdfUniform);

		vector<vector<float>> vpdfSkins;
		for(i=0; i<vbuckets.size(); ++i)
		{
			vpdfSkins.push_back(computePDF(vbuckets[i], minVal, maxVal, numBins));
		}
		float avJSD = 0;
		for(i=0; i<vpdfSkins.size(); ++i)
		{
			float jsd = computeJSD(vpdfSkins[i], pdfUniform);
			avJSD += jsd;
		}
		avJSD /= (float)vpdfSkins.size();
		features[0].Value = jsdAll - avJSD;
		features[1].Value = jsdAll;
	}

	return features;
}
