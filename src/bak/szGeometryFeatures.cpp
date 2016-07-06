
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
#include <szDistributionsTemplate.h>
#include <szDistributions.h>
#include <szEigen.h>
#include <szDistanceTransform.h>
//#include <szConvolutionTemplate.h>
#include <szMyNeighborOp.h>

inline int
NumGeometryFeatures()
{
	const int NumFeatures = 7;
	return NumFeatures;
}

vector<int>
removeSmallComponents(vector<int>& C, 
					  int numComponents,
					  int minCnt, 
					  int ndim, 
					  const int* dims)
{
	vector<int> vCnt(numComponents+1, 0);
	int i;
	for (i=0; i<C.size(); ++i)
	{
		if(C[i] <= numComponents)
		{
			vCnt[C[i]]++;
		}
	}
	//derive labels after removals of small components
	int label = 1;
	vector<int> vLabels(numComponents+1, 0);
	vector<int> vCnt2;
	for(i=1; i<vCnt.size(); ++i)
	{
		if(vCnt[i]>=minCnt)
		{
			vLabels[i] = label;
			label++;
			vCnt2.push_back(vCnt[i]);
		}
	}
	//re-label the components
	for (i=0; i<C.size(); ++i)
	{
		if(C[i]>0 && C[i] <= numComponents)
		{
			C[i] = vLabels[C[i]];
		}
	}
	return vCnt2; //this corresponds to the updated number of components
}

vector<CFeature>
computeGeometryFeatures(const vector<unsigned char>& B, //nodule segmentation
						const vector<unsigned char>& B0, //free trace segmentation
						const vector<unsigned char>& L, //foreground segmentation
						const vector<float>& D,         //distance map of B
						const vector<float>& D0,         //distance map of B0
						const vector<float>& Dl,         //distance map of L
						int xseed, int yseed, int zseed, //seed points
						int ndim,
						const int* dims) 
{
	int i, j, k, m;
	int nvoxels=numberOfElements(ndim,dims);

	vector<double> centroid = computeCentroid(B, ndim, dims);

	vector<float> vRad; //radii
	vector<float> vRough; //roughness
	for(k=0; k<dims[2]; ++k) 
	{
		for(i=0; i<dims[1]; ++i) 
		{
			for(j=0; j<dims[0]; ++j) 
			{
				if(GetData3(B, j, i, k, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					int surfaceCnt = 0;
					double sx = 0, sy = 0, sz = 0; //compute the surface normal at boundary 
					for(m=0; m<6; ++m)
					{
						int x2=j+XOffset6[m];
						int y2=i+YOffset6[m];
						int z2=k+ZOffset6[m];
						if(!GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)0))
						{
							sx += XOffset6[m];
							sy += YOffset6[m];
							sz += ZOffset6[m];
							surfaceCnt ++;
						}
					}
					if(surfaceCnt)
					{
						double dx = (j-centroid[0]);
						double dy = (i-centroid[1]);
						double dz = (k-centroid[2]);
						double length1 = sqrt(dx*dx + dy*dy + dz*dz);
						double length2 = sqrt(sx*sx + sy*sy + sz*sz);//!!!!!
						double ang = 0;
						if (length1>0 && length2>0)
							vRough.push_back((dx*sx + dy*sy + dz*sz) / (length1*length2));//!!!!!
						vRad.push_back(length1);
					}
				}
			}
		}
	}
	vector<double> vMomentRad = computeMoments(vRad, 2);
	vector<double> vMomentRough = computeMoments(vRough, 2);

	//Sphericity features
	float maxSphericity = MaximumSphericityValue(2, 2, 2);
	float spval=EvaluateSphericityValue(D, xseed, yseed, zseed, 2,2,2,dims) / maxSphericity;
	float spval0=EvaluateSphericityValue(D0, xseed, yseed, zseed, 2,2,2,dims) / maxSphericity;
	float spvalL=EvaluateSphericityValue(Dl, xseed, yseed, zseed, 2,2,2,dims) / maxSphericity;

	vector<CFeature> features(NumGeometryFeatures());
	features[0].Value = spval;
	features[0].Name = "SphericityD";
	features[1].Value = spval0;
	features[1].Name = "SphericityD0";
	features[2].Value = spvalL;
	features[2].Name = "SphericityDL";
	features[3].Value = spval* spval0 * spvalL;
	features[3].Name = "SphericityProd";
	if(vMomentRad[0])
		features[4].Value = 1.0 - sqrt(vMomentRad[1]) / vMomentRad[0]; //!!!!!
	features[4].Name = "SegRadiiSpread";
	double PI = 4*atan(1.0);
	features[5].Value = vMomentRough[0] - vMomentRough[1];
	features[5].Name = "SegSurfaceRoughness";
	features[6].Value = features[4].Value * features[5].Value;
	features[6].Name = "SegRadiiRoughnessProd";

	return features;
}

vector<CFeature>
computeGeometryFeaturesNew(const vector<unsigned char>& B, //nodule segmentation
						   const vector<unsigned char>& B0, //free trace result
						   const vector<unsigned char>& L, //foreground including the segmentation
						   const vector<float>& D,
						   const vector<float>& D0,
						   const vector<float>& Dl,
						   int xseed, int yseed, int zseed,
						   const int* dims) 
{
	int ndim = 3;
	vector<CFeature> features;
	int nvoxels = numberOfElements(ndim, dims);

	int numBins = 20;
	{
		int i;
		float maxVal = 1;
		int x0 = xseed, y0 = yseed, z0 = zseed;
		for(i=0; i<D.size(); ++i)
		{
			if(maxVal < D[i])
			{
				maxVal = D[i];
				Ind2Sub(x0, y0, z0, i, dims);
			}
		}
		float dApart = 1;
		for(i=0; i<D.size(); ++i)
		{
			if(B[i])
			{
				int x1, y1, z1;
				Ind2Sub(x1, y1, z1, i, dims);
				float d = sqrt((float)(x1-x0)*(x1-x0)+(y1-y0)*(y1-y0)+(z1-z0)*(z1-z0));
				if(d > dApart)
					dApart = d;
			}
		}
		double daRatio = maxVal / dApart;
		CFeature feature;
		feature.Name = "DistanceApartRatio";
		feature.Value = daRatio;
		features.push_back(feature);
	}
	{   
		int i, j, k;
		float maxSp = 0;
		float maxSp0 = 0;
		int sphNbr[3] = {1, 1, 1};
		float maxSphericity = MaximumSphericityValue(1,1,1);
		for(i=0; i<dims[2]; ++i) 
		{
			for(j=0; j<dims[1]; ++j) 
			{
				for(k=0; k<dims[0]; ++k) 
				{
					if((GetData3(B0, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0) ||
						GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0)))
					{
						float spval = EvaluateSphericityValue(D0,k,j,i,sphNbr[0], sphNbr[1], sphNbr[2],dims)/maxSphericity;
						if(GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
						{
							if(spval > maxSp)
							{
								maxSp = spval;
							}
						}
						else
						{
							if(spval > maxSp0)
							{
								maxSp0 = spval;
							}
						}
					}
				}
			}
		}
		CFeature feature1("MaxSegSphericity", maxSp);
		CFeature feature2("SegSphericityRatio", 0);
		if(maxSp>0)
			feature2.Value = 1 - maxSp0/maxSp;

		features.push_back(feature1);
		features.push_back(feature2);
	}
	{   
		//compute KSD between D(segmentation) and D(Foreground @segmentation).
		//the two should have similar distributions.
		int i;
		vector<float> Db(Dl.size(), 0); //D(Foreground @segmentation)
		for(i=0; i<L.size(); ++i)
		{
			if(B[i])
				Db[i] = Dl[i];
		}
		//find the maximum value of in the two distance maps 
		float maxVal = 1;
		for(i=0; i<Db.size(); ++i)
		{
			maxVal = Max(maxVal, D[i]);
			maxVal = Max(maxVal, Db[i]);
		}
		vector<float> vedges;
		for(i=1; i<=(int)maxVal; ++i)
		{
			vedges.push_back((float)i);
		}
		vedges.push_back((float)i);

		vector<float> cb = computeCDF(Db, vedges);
		vector<float> cc = computeCDF(D, vedges);
		double ksd = 1.0 - computeKSD(cb, cc);
		CFeature feature("DmapDisparity", ksd);
		features.push_back(feature);
	}
	return features;
}

vector<CFeature>
computeErosionFeatures(const vector<unsigned char>& B, //segmentation
					   const vector<unsigned char>& B0, //free trace segmentation
					   const vector<float>& D0,			//distance map of free trace segmentation
					   const int* dims)
{
	vector<CFeature> features;
					
	int i, j, k;
	vector<unsigned char> J(B.size(), 0); //joint area.
	int nj=0;
	int ns=0;
	for(i=0; i<dims[2]; ++i) 
	{
		for(j=0; j<dims[1]; ++j) 
		{
			for(k=0; k<dims[0]; ++k) 
			{
				if(GetData3(B0, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0) && 
					!GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					for(int m=0; m<NumNeighbors6; ++m)
					{
						int i2 = i + ZOffset6[m];
						int j2 = j + YOffset6[m];
						int k2 = k + XOffset6[m];
						if(GetData3(B, k2, j2, i2, dims[0], dims[1], dims[2], (unsigned char)0))
						{
							SetData3(J, k, j, i, dims[0], dims[1], dims[2], (unsigned char)1);
							nj++;
						}
					}
				}
				if(onSurface3(B, k, j, i, dims))
				{
					ns++;
				}
			}
		}
	}

	float dmaxJ = 0, dmaxB = 0;
	for(i=0; i<dims[2]; ++i) 
	{
		for(j=0; j<dims[1]; ++j) 
		{
			for(k=0; k<dims[0]; ++k) 
			{
				if(GetData3(J, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					float dval = GetData3(D0, k, j, i, dims[0], dims[1], dims[2], (float)0);
					dmaxJ = Max(dval, dmaxJ);
				}
				if(GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					float dval = GetData3(D0, k, j, i, dims[0], dims[1], dims[2], (float)0);
					dmaxB = Max(dval, dmaxB);
				}
			}
		}
	}
	CFeature feature1("JointBoundaryFitness1");
	CFeature feature2("JointBoundaryFitness2");
	CFeature feature3("JointBoundaryFitness3");
	CFeature feature4("JointBoundaryFitness4");

	if(dmaxB>0)
	{
		feature1.Value = exp(-dmaxJ/(dmaxB*dmaxB));
		feature2.Value = 1.0/(1.0 + exp(dmaxJ-dmaxB));
		feature3.Value = 1.0/(1.0 + exp(dmaxJ*dmaxJ-dmaxB*dmaxB));
		if(ns)
			feature4.Value = 1.0 - (double)nj/(double)ns;
	}
	features.push_back(feature1);
	features.push_back(feature2);
	features.push_back(feature3);
	features.push_back(feature4);
	return features;
}
