#ifdef MEX_DLL
#include <mex.h>
#endif
#include <szMexUtilityTemplate.h>
#include <szDistributions.h>
#include <szFeatureClass.h>
#include <szPBSFeatures.h>
#include <szDistributionsTemplate.h>

template<class Item>
vector<int>
computeFrequency(const vector<Item>& values,
				 const Item& minVal,
				 const Item& maxVal,
				 int numBins)
{
	vector<int> count(numBins, 0);
	if(values.empty())
	{
		return count;
	}

	int i;
	if(minVal == maxVal)
	{
		//obtain minVal and maxVal from the data.
		Item minVal = values[0];
		Item maxVal = values[0];
		for(i=1; i<values.size(); ++i)
		{
			minVal = Min(minVal, values[i]);
			maxVal = Max(maxVal, values[i]);
		}
	}
	if(maxVal > minVal)
	{
		double interval = (maxVal - minVal) / numBins;
		for(i=0; i<values.size(); ++i)
		{
			int k = (values[i] - minVal) / interval;
			if(k>=0 && k<numBins)
			{
				count[k]++;
			}
		}	
	}
	else
	{
		count[0] = values.size();
	}

	return count;
}

template<class Item>
vector<int>
computeFrequency(const vector<Item>& values,
				 const vector<Item>& vedges)
{
	if(vedges.size()<2) //needs at least two bin edges
	{
		return vector<int>(0);
	}
	vector<int> vcount(vedges.size() - 1, 0);
	if(values.empty())
	{
		return vcount;
	}

	int i, j;
	for(i=0; i<values.size(); ++i)
	{
		if(vedges[0] > values[i])
			continue;

		int id = binarySearch(vedges, values[i]);
		if(values[i] < vedges[id] && id>0)
		{
			id--;
		}
		vcount[id]++;
	}	

	return vcount;
}

template<class Item>
vector<float>
computePDF(const vector<Item>& values,
		   const Item& minVal,
		   const Item& maxVal,
		   int numBins)
{
	if(values.empty())
	{
		return vector<float>(numBins, 0);
	}
	vector<int> vhist = computeFrequency(values, minVal, maxVal, numBins);
	vector<float> pdf = computePDF(vhist);

	return pdf;
}

template<class Item>
vector<float>
computePDF(const vector<Item>& values,
		   const vector<Item>& vedges)
{
	vector<int> vhist = computeFrequency(values, vedges);
	vector<float> pdf = computePDF(vhist);

	return pdf;
}

template<class Item>
vector<float>
computeCDF(const vector<Item>& values,
		   const Item& minVal,
		   const Item& maxVal,
		   int numBins)
{
	vector<float> pdf = computePDF(values, minVal, maxVal, numBins);
	vector<float> cdf = computeCDF(pdf);

	return cdf;
}

template<class Item>
vector<float>
computeCDF(const vector<Item>& values,
		   const vector<Item>& vedges)
{
	vector<float> pdf = computePDF(values, vedges);
	vector<float> cdf = computeCDF(pdf);

	return cdf;
}

template<class Item>
double
JSD_Features(const vector<Item>& A,
			 const vector<Item>& B,
			 int numBins)
{
	int i;
	Item minVal = allvals[0];
	Item maxVal = allvals[0];
	for(i=1; i<A.size(); ++i)
	{
		minVal = Min(minVal, A[i]);
		maxVal = Max(maxVal, A[i]);
	}
	for(i=1; i<B.size(); ++i)
	{
		minVal = Min(minVal, B[i]);
		maxVal = Max(maxVal, B[i]);
	}

	vector<float> pdfA = computePDF(A, minVal, maxVal, numBins);
	vector<float> pdfB = computePDF(A, minVal, maxVal, numBins);

	return computeJSD(pdfA, pdfB);
}

template<class Item1, class Item2>
float 
computeCorrelation(const vector<Item1>& va,
				   const vector<Item2>& vb,
				   int length)
{
	if(length <= 0)
	{
		length = Min(va.size(), vb.size());
	}
	if(length == 0)
	{
		return (float)0;
	}

	int i;
	double meanA=0, meanB=0;
	double engA=0, engB=0;
	for(i=0; i<length; ++i)
	{
		meanA += va[i];
		meanB += vb[i];
		engA += va[i]*va[i];
		engB += vb[i]*vb[i];
	}
	meanA /= length;
	meanB/= length;
	double stdA = sqrt(engA/length - meanA*meanA);
	double stdB = sqrt(engB/length - meanB*meanB);
	double corr = 0;
	if(stdA>0 && stdB>0)
	{
		for(i=0; i<length; ++i)
		{
			corr += (va[i] - meanA)*(vb[i] - meanB);
		}
		corr /= (stdA*stdB*length);
	}
	return (float)corr;
}
