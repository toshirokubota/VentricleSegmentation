
#include <string>
#include <cmath>
#include <algorithm>
#include <cstdio>
using namespace std;

#include <szMiscOperations.h>
#include <szMexUtility.h>	
#include <szMexUtilityTemplate.h>
#include <szDistanceTransform.h>
#include <szDistributionsTemplate.h>
#include <szPBSFeatures.h>

inline int
NumPBSFeatures()
{
	const int NumFeatures = 6;
	return NumFeatures;
}

bool
PBS_IndexVector(vector<int>& S, 
				const vector<float>& Dsq, 
				int nbins)
{
	float dmax;
	getMaximum(Dsq, dmax);
	if(dmax<=0) //no segmentation is found
		return false;

	vector<float> vedges((int)dmax+1, 0);
	int i;
	for(i=0; i<vedges.size(); ++i)
	{
		vedges[i] = (float)i+1;
	}
	vector<float> vcdf = computeCDF(Dsq, vedges);
	vector<int> vlabel(vcdf.size());
	float delta = 1.0 / nbins;
	float thres = delta;
	int lb = 1;
	for(i=0; i<vlabel.size(); ++i)
	{
		if(vcdf[i] <= thres)
		{
			vlabel[i] = lb;
		}
		else
		{
			if(vcdf[i]-thres>delta/2)
				vlabel[i] = lb+1;
			else
				vlabel[i] = lb;

			lb++;
			thres += delta;
		}
	}
	for(i=0; i<Dsq.size(); ++i)
	{
		if(Dsq[i])
		{
			S[i] = vlabel[(int)Dsq[i] - 1];
		}
		else
		{
			S[i] = 0;
		}
	}

	return true;
}

bool
PBS_IndexVectorWrapperShell(vector<int>& S, 
							const vector<float>& D, 
							int nshells)
{
	int i;
	vector<float> D2(D.size());
	for(i=0; i<D.size(); ++i)
	{
		D2[i] = D[i] * D[i];
	}
	bool bRet = PBS_IndexVector(S, D2, nshells);
	return bRet;
}

bool
PBS_IndexVectorWrapperWedge(vector<int>& S, //index volume
							const vector<unsigned char>& B, //segmentation
							double x, double y, double z, //segmentation center
							int nwedges,
							int nbins,
							const int* dims)
{
	vector<float> G(B.size());
	int i, j, k;
	double PI = 3.14159265;
	int count = 0;
	for(i=0; i<dims[2]; ++i)
	{
		double dz = i-z;
		for(j=0; j<dims[1]; ++j)
		{
			double dy = j-y;
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					double dx = k-x;
					double lenXY = sqrt(dx*dx + dy*dy);
					double phi = atan2(dy, dx);
					double psi = atan2(dz, lenXY);
					int nphi = Round((phi + PI)/nwedges);
					int npsi = Round((psi + PI)/nwedges);
					SetData3(S, k, j, i, dims[0], dims[1], dims[2], (int)(npsi * nwedges + nphi + 1));
					count++;
				}
				else
				{
					SetData3(S, k, j, i, dims[0], dims[1], dims[2], (int)0);
				}
			}
		}
	}
	if(count)
		return true;
	else
		return false;
}

/*vector<CFeature>
computePBSFeatures(const vector<unsigned short>& A, 
				   const vector<unsigned char>& B, 
				   const vector<float>& Ds, 
				   int x, int y, int z,
				   const int* dims)
{
	vector<CFeature> vFeatures;

	//PBS features
	int i;
	vector<float> D2(Ds.size());
	for(i=0; i<Ds.size(); ++i)
	{
		D2[i] = Ds[i] * Ds[i];
	}
	//vector<int> S(D.size(), 0);
	if(numMoreThan(B, (unsigned char) 0)) //at least one voxel is required.
	{
		vector<int> S(A.size(), 0); 
		int nshells = 20;
		PBS_IndexVector(S, D2, nshells);

		int nbins = 20;
		vector<CFeature> PBS = PBS_ShellFeatures(A, S, nbins, "Gray", dims);
		vFeatures.insert(vFeatures.end(), PBS.begin(), PBS.end());
		printf("PBS: %f, %f\n", PBS[0], PBS[1]);

		vector<float> G(A.size(), 0);
		gradientMagnitude(G, A, B, dims);
		vector<CFeature> PBSG = PBS_ShellFeatures(G, S, nbins, "GradMag", dims);
		vFeatures.insert(vFeatures.end(), PBSG.begin(), PBSG.end());
		printf("PBSG: %f, %f\n", PBSG[0], PBSG[1]);

		vector<float> T(A.size(), 0);
		gradientAngle(T, A, B, x, y, z, dims);
		vector<CFeature> PBSA = PBS_ShellFeatures(T, S, nbins, "GradAng", dims);
		vFeatures.insert(vFeatures.end(), PBSA.begin(), PBSA.end());
		printf("PBSA: %f, %f\n", PBSA[0], PBSA[1]);
	}
	else
	{
		for(int k=0; k<6; ++k)
		{
			vFeatures.push_back(0);
		}
	}

	return vFeatures;
}*/