inline
vector<CFeature>
correlationFeatures(const char* name)
{	
	int numFeatures = 2;
	vector<CFeature> features(numFeatures);
	features[0].Name = GetCorrelationFeatureName1(name);
	features[1].Name = GetCorrelationFeatureName2(name);

	return features;
}

template<class Item1, class Item2>
vector<CFeature>
computeCorrelationFeatures(const vector<Item1>& A, //image 1
						   const vector<Item2>& B,  //image 2
						   const vector<unsigned char>& M, //mask 
						   const char* name, 
						   const int* dims)
{
	vector<CFeature> features = correlationFeatures(name);
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

	features[0].Value = innerP;
	features[1].Value = corr;

	return features;
}
