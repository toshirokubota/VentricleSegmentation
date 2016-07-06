#include <mex.h>
#include <szMexUtilityTemplate.h>
#include <szDistributions.h>

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

/*
for each voxel in B, find the bin number according to vedges and 
store the number in S.
*/
template<class Item>
bool
IndexVector(vector<int>& S, 
			const vector<Item>& B, 
			const vector<Item>& vedges)
{
	int i;
	for(i=0; B.size(); ++i)
	{
		int id = binarySearch(vedges, B[i]);
		S[i] = id;
	}

	return true;
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
	//printf("muA = %f, muB = %f, stdA = %f, stdB = %f\n", 
	//	meanA, meanB, stdA, stdB);
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

template<class Item>
vector<CFeature>
PBS_ShellFeatures(const vector<Item>& A,
				  const vector<int>& D,
				  int numBins,
				  const char* name,
				  const int* dims)
{
	vector<CFeature> features;
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
	string name1 = GetPBSFeatureName1(name);
	CFeature feature1(name1);
	string name2 = GetPBSFeatureName2(name);
	CFeature feature2(name2);
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
		//printf("PBS_ShellFeatures: #skins = %d, avJSD = %f\n", vpdfSkins.size(), avJSD);
		feature1.Value = jsdAll - avJSD;
		feature2.Value = jsdAll;
	}
	features.push_back(feature1);
	features.push_back(feature2);

	return features;
}

template<class Item>
void
gradientMagnitude(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims)
{
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					Item c = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (Item)0);
					Item w = GetData3(A, k-1, j, i, dims[0], dims[1], dims[2], c);
					Item e = GetData3(A, k+1, j, i, dims[0], dims[1], dims[2], c);
					Item n = GetData3(A, k, j-1, i, dims[0], dims[1], dims[2], c);
					Item s = GetData3(A, k, j+1, i, dims[0], dims[1], dims[2], c);
					Item f = GetData3(A, k, j, i-1, dims[0], dims[1], dims[2], c);
					Item b = GetData3(A, k, j, i+1, dims[0], dims[1], dims[2], c);
					float gr = sqrt((float)(w-e)*(w-e) + (n-s)*(n-s) + (f-b)*(f-b));
					SetData3(G, k, j, i, dims[0], dims[1], dims[2], gr);
				}
			}
		}
	}
	return;
}

template<class Item>
void
laplacianResponse(vector<float>& G,
				  const vector<Item>& A,
				  const vector<unsigned char>& S,
				  const int* dims)
{
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					Item c = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (Item)0);
					Item w = GetData3(A, k-1, j, i, dims[0], dims[1], dims[2], c);
					Item e = GetData3(A, k+1, j, i, dims[0], dims[1], dims[2], c);
					Item n = GetData3(A, k, j-1, i, dims[0], dims[1], dims[2], c);
					Item s = GetData3(A, k, j+1, i, dims[0], dims[1], dims[2], c);
					Item f = GetData3(A, k, j, i-1, dims[0], dims[1], dims[2], c);
					Item b = GetData3(A, k, j, i+1, dims[0], dims[1], dims[2], c);
					float lp = (float) c - ((float)w + e + s + n + f + b) / 6.0;
					SetData3(G, k, j, i, dims[0], dims[1], dims[2], lp);
				}
			}
		}
	}
	return;
}

template<class Item>
void
gradientAngle(vector<float>& G,
			  const vector<Item>& A,
			  const vector<unsigned char>& S,
			  int x, int y, int z,
			  const int* dims)
{
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		int dz = i-z;
		for(j=0; j<dims[1]; ++j)
		{
			int dy = j-y;
			for(k=0; k<dims[0]; ++k)
			{
				int dx = k-x;
				if(GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					Item c = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (Item)0);
					Item w = GetData3(A, k-1, j, i, dims[0], dims[1], dims[2], c);
					Item e = GetData3(A, k+1, j, i, dims[0], dims[1], dims[2], c);
					Item n = GetData3(A, k, j-1, i, dims[0], dims[1], dims[2], c);
					Item s = GetData3(A, k, j+1, i, dims[0], dims[1], dims[2], c);
					Item f = GetData3(A, k, j, i-1, dims[0], dims[1], dims[2], c);
					Item b = GetData3(A, k, j, i+1, dims[0], dims[1], dims[2], c);
					float ang = 0;
					if(dz == 0 && dy == 0 && dx==0)
					{
						ang = 0;
					}
					else if(w==e && s==n && f==b)
					{
						ang = 0;
					}
					else
					{
						float l1 = sqrt((float)dx*dx + dy*dy + dz*dz);
						float l2 = sqrt((float)(w-e)*(w-e) + (s-n)*(s-n) + (f-b)*(f-b));
						ang = (dx*(e-w)+dy*(s-n)+dz*(b-f))/(l1*l2);
					}
					SetData3(G, k, j, i, dims[0], dims[1], dims[2], ang);
				}
			}
		}
	}
	return;
}

template<class Item1, class Item2>
vector<CFeature>
computeCorrelationFeatures(const vector<Item1>& A, //image 1
						   const vector<Item2>& B,  //image 2
						   const vector<unsigned char>& M, //mask 
						   const char* name, 
						   const int* dims)
{
	vector<CFeature> features;
	int i, j, k;
	//inner product feature
	double length1 = 0;
	double length2 = 0;
	double dotp = 0;
	int numSegVol = 0;
	for(i=0; i<A.size(); ++i)
	{
		if(M[i])
		{
			length1 += A[i] * A[i];
			length2 += B[i] * B[i];
			dotp += A[i] * B[i];
			numSegVol++;
		}
	}
	double innerP = 0;
	double corr = 0;
	if(numSegVol)
	{
		if(length1>0 && length2>0)
		{
			innerP = dotp /(sqrt(length1)*(sqrt(length2)));
		}
		//correlation feature
		vector<float> valsA(numSegVol);
		vector<float> valsB(numSegVol);
		int count = 0;
		for(i=0; i<A.size(); ++i)
		{
			if(M[i])
			{
				valsA[count] = A[i];
				valsB[count] = B[i];
				count++;
			}
		}
		corr = computeCorrelation(valsA, valsB);
	}

	string name1 = GetCorrelationFeatureName1(name);
	CFeature feature1(name1, innerP);
	string name2 = GetCorrelationFeatureName2(name);
	CFeature feature2(name2, corr);

	features.push_back(feature1);
	features.push_back(feature2);

	return features;
}