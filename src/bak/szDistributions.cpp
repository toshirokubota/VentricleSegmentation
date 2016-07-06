#include <szDistributions.h>
#include <szMexUtility.h>
#include <cmath>
using namespace std;

vector<float>
computePDF(const vector<int>& vhist)
{
	int numBins = vhist.size();
	vector<float> pdf(numBins, 0);

	int i;
	int count = 0;
	for(i=0; i<numBins; ++i)
	{
		count += vhist[i];
	}
	//when this happen, we set PDF to all zero.
	//Although, this violates the probability axiom, it appears to be the most
	//harmless way to take care of such ill-condition. 
	if(count) 
	{
		for(i=0; i<numBins; ++i)
		{
			pdf[i] = (float)vhist[i] / count;
		}
	}
	return pdf;
}

vector<float>
computeCDF(const vector<float>& pdf)
{
	vector<float> cdf(pdf.size(), 0);
	int i;
	float cnt = 0;
	for(i=0; i<cdf.size(); ++i)
	{
		cdf[i] = cnt + pdf[i];
		cnt += pdf[i];
	}

	return cdf;
}

vector<float>
computeCDF(const vector<int>& vhist)
{
	vector<float> pdf = computePDF(vhist);
	vector<float> cdf = computeCDF(pdf);

	return cdf;
}

float 
computeKSD(const vector<float>& ca,
		   const vector<float>& cb)
{
	int nbins = ca.size(); 
	if(nbins == 0)
		return 0;

	float ksd = 0;
	int i;
	for(i=0; i<nbins; ++i)
	{
		ksd = Max(ksd, Abs(ca[i]-cb[i]));
	}

	return ksd;
}

float
computeJSD(const vector<float>& pa,
		   const vector<float>& pb)
{
	int nbins = pa.size(); 
	if(nbins == 0)
		return 0;

	float jsd = 0;
	int i;
	for(i=0; i<nbins; ++i)
	{
		float la, lb;
		if(pa[i] <= 0 && pb[i] <= 0)
		{
			continue;
		}
		if(pa[i] <= 0)
		{
			la = 0;  
		}
		else
		{
			la = log(pa[i]);
		}
		if(pb[i] <= 0)
		{
			lb = 0;  
		}
		else
		{
			lb = log(pb[i]);
		}
		jsd -= pa[i]*la + pb[i]*lb - (pa[i]+pb[i])*log(pa[i]+pb[i]);
	}

	return jsd/2.0;
}
