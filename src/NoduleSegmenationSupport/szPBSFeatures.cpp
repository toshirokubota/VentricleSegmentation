
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

string
GetPBSFeatureName1(const char* name0)
{
	string name = "PBSDif";
	name += name0;
	return name;
}

string
GetPBSFeatureName2(const char* name0)
{
	string name = "PBS";
	name += name0;
	return name;
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
