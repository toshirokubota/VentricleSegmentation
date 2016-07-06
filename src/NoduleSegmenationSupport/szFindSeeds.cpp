
#ifdef MEX_DLL
#include <mex.h>
#endif
#include <string.h>
#include <math.h>
#include <algorithm>
using namespace std;
#include <stdio.h>
#include <szDefaultParam.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szConnectedComponent.h>
#include <szDistanceTransform.h>
#include <szLocalExtrema.h>
#include <szMiscOperations.h>
#include <szParticle.h>


/*
Another version that is very similar to the other FindSeedsThres.
The approach is exactly the same and only the interface is different.
Instead of returning a vector of seed points, it stores the seeds in\
vseeds and returns the maximum spherical value.
*/
vector<int>
FindSeedsThres(const vector<unsigned char>& L,  //INPUT - foreground segmentation
			   vector<float>& D,          //OUTPUT - distance map
			   vector<float>& Sp,          //OUTPUT - sphericity map
			   vector<unsigned char>& P,  //OUTPUT - local maxima of the sphericity map
			   vector<float>& vSpherical, //OUTPUT - spherical values of seeds
			   int w,
			   int wcore, 
			   float selection_threshold,
			   const int* dims) 
{

	int i, j, k;
	int ndim=3; //3D is assumed
	DistanceTransformFInv(D, L, DistanceTransformMode_Euclid, ndim, dims);
	//printf("Done with Distance Transform.\n");

	int x0=dims[0]/2, y0=dims[1]/2, z0=dims[2]/2;
	int sphNbr0[3] = {1, 1, 1};
	float maxSphericity1 = MaximumSphericityValue(sphNbr0[0],sphNbr0[1],sphNbr0[2]);
	for(i=0; i<dims[2]; ++i) 
	{ 
		for(j=0; j<dims[1]; ++j) 
		{ 
			for(k=0; k<dims[0]; ++k) 
			{ 
				float val = EvaluateSphericityValue(D,k,j,i,sphNbr0[0],sphNbr0[1],sphNbr0[2],dims)/maxSphericity1;
				SetData3(Sp,k,j,i,dims[0],dims[1],dims[2],val);
			}
		}
	}
	vector<int> nbh = MakeEightNeighborhood(ndim,dims);
	LocalMaximum(P,Sp,L,nbh,false,ndim,dims);
	int sphNbr[3] = {2, 2, 2};
	//printf("Done with Local Maximum.\n");

	//mark local maximum within a confined window located at the center
	vector<float> vvalues;
	vector<int> vind;
	float maxSphericity = MaximumSphericityValue(sphNbr[0], sphNbr[1], sphNbr[2]);
	//printf("Maximum possible sphericity is %f\n", maxSphericity);
	for(i=z0-w; i<=z0+w; ++i) 
	{ 
		for(j=y0-w; j<=y0+w; ++j) 
		{ 
			for(k=x0-w; k<=x0+w; ++k) 
			{ 
				if(GetData3(P,k,j,i,dims[0],dims[1],dims[2],(unsigned char) 0)) 
				{
					float dval = GetData3(D,k,j,i,dims[0],dims[1],dims[2],(float) 0);
					float spval = GetData3(Sp,k,j,i,dims[0],dims[1],dims[2],(float) 0);
					//if(dval>1 && spval>0.0) //d>1 == cannot be at the boarder of the foreground
					if(spval>0.0) //!!!TK 02/26/2008
						//if(!onSurfaceOnSlice(L,i,j,k,dims) && spval>0.0) //d>1 == cannot be at the boarder of the foreground
					{
						float val=EvaluateSphericityValue(D,k,j,i,sphNbr[0], sphNbr[1], sphNbr[2],dims);
						//printf("(%d, %d, %d): %f, %f, %f\n", k, j, i, dval, spval, val);
						if(val>0)  //valid value should be positive
						{
							/*float dChecker = Max(Abs(i-z0),Max(Abs(j-y0),Abs(k-x0)));
							if(dChecker>wcore)  //weight the value if it is too far away from the center
							{
							val*=(1.-(dChecker-wcore)/wcore);
							}*/
							val /= maxSphericity;
							vvalues.push_back(val);
							vind.push_back(Sub2Ind(k,j,i,dims));
							//SetData3(Sp,k,j,i,dims[0],dims[1],dims[2],(float)(val + 1.0));
							//printf("Sphericity: (%d, %d, %d) %f\n", k,j,i,val);
						}
					}
				}
			}
		}
	}
	vector<int> vsub;
	if(vvalues.size()) 
	{
		vector<int> vindex = indexedSort(vvalues);
		selection_threshold = Min(vvalues[vvalues.size()-1]*0.95, selection_threshold);
		//printf("max sphericity: %f\n", vvalues[vvalues.size()-1]); //!!!!!
		if(vvalues[vvalues.size()-1] >= DefaultMinimumValidSphericity)
		{
			for(i=0; i<vvalues.size() && (i==0 || vvalues[vvalues.size()-i-1] >= selection_threshold); ++i) 
			{
				int x,y,z;
				Ind2Sub(x,y,z,vind[*(vindex.end()-1-i)],dims);
				vsub.push_back(x);
				vsub.push_back(y);
				vsub.push_back(z);
				vSpherical.push_back(vvalues[vvalues.size()-i-1]);
				printf("Seed at (%d %d %d) %2.2f\n", x+1, y+1, z+1, vvalues[vvalues.size()-i-1]);
			}
		}
	}
	return vsub;
}

/*
Let a region grow from each seed point for a specified number of iterations (niter).
The map provides regions in which seed points are to be sought.
M stores the regions (as non-zero value).
*/
void
GrowRegionFromClickPoints(vector<unsigned char>& M,
						  const vector<unsigned char>& L,
						  vector<CParticle> vclick, 
						  int niter, 
						  const int* dims)
{
	for(int n=0; n<niter; ++n)
	{
		vector<CParticle> vnext;
		for(int c=0; c<vclick.size(); ++c)
		{
			int x0 = vclick[c].m_X;
			int y0 = vclick[c].m_Y;
			int z0 = vclick[c].m_Z;
			for(int i=z0-1; i<=z0+1; ++i) 
			{ 
				for(int j=y0-1; j<=y0+1; ++j) 
				{ 
					for(int k=x0-1; k<=x0+1; ++k) 
					{
						if(GetData3(L,k,j,i,dims[0],dims[1],dims[2],(unsigned char) 0)) 
						{
							if(!GetData3(M,k,j,i,dims[0],dims[1],dims[2],(unsigned char) 0)) 
							{
								vnext.push_back(CParticle(k, j, i));
								SetData3(M,k,j,i,dims[0],dims[1],dims[2],(unsigned char) 1);
							}
						}
					}
				}
			}
		}
		vclick = vnext;
	}
}

/*
This uses user supplied click points.  
*/
vector<CParticle>
FindSeedsThres(const vector<unsigned char>& L,  //INPUT - foreground segmentation
			   vector<float>& D,          //OUTPUT - distance map
			   vector<float>& Sp,          //OUTPUT - sphericity map
			   vector<unsigned char>& P,  //OUTPUT - local maxima of the sphericity map
			   vector<float>& vSpherical, //OUTPUT - spherical values of seeds
			   const vector<CParticle>& vclick, //INPUT - click points
			   int w,
			   float selection_threshold,
			   const int* dims) 
{

	int i, j, k;
	int ndim=3; //3D is assumed
	DistanceTransformFInv(D, L, DistanceTransformMode_Euclid, ndim, dims);
	int sphNbr0[3] = {1, 1, 1};
	float maxSphericity1 = MaximumSphericityValue(sphNbr0[0],sphNbr0[1],sphNbr0[2]);
	vector<unsigned char> M(L.size(), 0);
	GrowRegionFromClickPoints(M, L, vclick, w, dims);

	for(i=0; i<dims[2]; ++i) 
	{ 
		for(j=0; j<dims[1]; ++j) 
		{ 
			for(k=0; k<dims[0]; ++k) 
			{ 
				float val = EvaluateSphericityValue(D,k,j,i,sphNbr0[0],sphNbr0[1],sphNbr0[2],dims)/maxSphericity1;
				SetData3(Sp,k,j,i,dims[0],dims[1],dims[2],val);
			}
		}
	}
	vector<int> nbh = MakeEightNeighborhood(ndim,dims);
	LocalMaximum(P,Sp,L,nbh,false,ndim,dims);
	int sphNbr[3] = {2, 2, 2};

	//mark local maximum within a confined window located at a click point
	vector<float> vvalues;
	vector<CParticle> vseedP;
	float maxSphericity = MaximumSphericityValue(sphNbr[0], sphNbr[1], sphNbr[2]);
	for(int c=0; c<vclick.size(); ++c)
	{
		int x0 = vclick[c].m_X;
		int y0 = vclick[c].m_Y;
		int z0 = vclick[c].m_Z;
		for(i=z0-w; i<=z0+w; ++i) 
		{ 
			for(j=y0-w; j<=y0+w; ++j) 
			{ 
				for(k=x0-w; k<=x0+w; ++k) 
				{ 
					if(GetData3(P,k,j,i,dims[0],dims[1],dims[2],(unsigned char) 0) &&
						GetData3(M,k,j,i,dims[0],dims[1],dims[2],(unsigned char) 0)) 
					{
						float dval = GetData3(D,k,j,i,dims[0],dims[1],dims[2],(float) 0);
						float spval = GetData3(Sp,k,j,i,dims[0],dims[1],dims[2],(float) 0);
						if(spval>0.0) //!!!TK 02/26/2008
						{
							CParticle p(k, j, i);
							if(find(vseedP.begin(), vseedP.end(), p) == vseedP.end())
							{
								float val=EvaluateSphericityValue(D,k,j,i,sphNbr[0], sphNbr[1], sphNbr[2],dims) / maxSphericity;
								if(val>0)  //valid value should be positive
								{
									vvalues.push_back(Max(val, spval)); //a small nodule can be penalized by using a bigger delta. 01/11/09
									vseedP.push_back(p);
									//printf("FindSeeds: candidate at (%d %d %d) %2.2f\n", p.m_X+1, p.m_Y+1, p.m_Z+1, val);
								}
							}
						}
					}
				}
			}
		}
	}
	vector<CParticle> vseedF;

	if(vvalues.size()) 
	{
		vector<int> vindex = indexedSort(vvalues);
		selection_threshold = Min(vvalues[vvalues.size()-1]*0.9, selection_threshold);
		printf("FindSeedsThres: the selection threshold is %f\n", selection_threshold);
		if(vvalues[vvalues.size()-1] >= DefaultMinimumValidSphericity)
		{
			for(i=0; i<vvalues.size() && (i==0 || vvalues[vvalues.size()-i-1] >= selection_threshold); ++i) 
			{
				CParticle p = vseedP[*(vindex.end()-1-i)];
				vseedF.push_back(p);
				vSpherical.push_back(vvalues[vvalues.size()-i-1]);
				printf("Seed at (%d %d %d) %2.2f\n", p.m_X+1, p.m_Y+1, p.m_Z+1, vvalues[vvalues.size()-i-1]);
			}
		}
	}
	return vseedF;
}

bool traceUphill(const vector<float>& D,
				 int x0, int y0, int z0,
				 int x, int y, int z,
				 int& xf, int& yf, int& zf,
				 float dthres,
				 const int* dims)
{
	if(x==0 || y==0 || z ==0 || x==dims[0]-1 || y==dims[1]-1 || z==dims[2]-1)
	{
		return false;
	}
	float dval1 = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	if(dval1 < 0)
	{
		return false;
	}
	float dmax = 0;
	vector<CParticle> pn;
	for(int i=0; i<NumNeighbors; ++i)
	{
		int x2 = x+XOffset[i];
		int y2 = y+YOffset[i];
		int z2 = z+ZOffset[i];
		float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
		if(dmax < dval2)
		{
			pn.clear();
			pn.push_back(CParticle(x2, y2, z2));
			dmax = dval2;
		}
		else if(dmax == dval2)
		{
			pn.push_back(CParticle(x2, y2, z2));
		}
	}
	//float dval2 = GetData3(D, xn, yn, zn, dims[0], dims[1], dims[2], (float)0);
	if(dmax > dval1)
	{
		for(int j=0; j<pn.size(); ++j)
		{
			int xn = pn[j].m_X;
			int yn = pn[j].m_Y;
			int zn = pn[j].m_Z;
			float dist = (float) Max(Abs(xn-x0), Max(Abs(yn-y0),Abs(zn-z0)));
			if(dist <= dthres)
			{
				if(traceUphill(D, x0, y0, z0, xn, yn, zn, xf, yf, zf, dthres, dims))
				{
					return true;
				}
			}
		}
	}
	else //the current position is on a local maximum
	{
		xf = x;
		yf = y;
		zf = z;
		//printf("traceHill: Found at %d, %d, %d (%f)\n", xf+1, yf+1, zf+1, dval1);
		return true;
	}
	return false;
}

vector<int>
refineSeeds(const vector<int>& vseeds,
			const vector<float>& D,
			const vector<unsigned char>& L,
			float dval0,
			const int* dims)
{
	int ndim = 3;
	int i;
	vector<int> refinedSeeds = vseeds;
	for(i=0; i<vseeds.size(); i+=3)
	{
		int x0 = vseeds[i];
		int y0 = vseeds[i+1];
		int z0 = vseeds[i+2];
		int xf, yf, zf;
		float dval;
		if(dval0 <= 0)
			dval = GetData3(D, x0, y0, z0, dims[0], dims[1], dims[2], (float)0);
		else
			dval = dval0;
		if(traceUphill(D, x0, y0, z0, x0, y0, z0, xf, yf, zf, ceil(dval), dims))
			//if(traceUphill(D, x0, y0, z0, x0, y0, z0, xf, yf, zf, 10, dims))
		{
			refinedSeeds[i] = xf;
			refinedSeeds[i+1] = yf;
			refinedSeeds[i+2] = zf;
			printf("Refine I: Seed %d is refined from (%d, %d, %d) to (%d, %d, %d)\n",
				i/3+1, vseeds[i]+1, vseeds[i+1]+1, vseeds[i+2]+1, xf+1, yf+1, zf+1);  
		}
	}
	return refinedSeeds;
}

/*
This version search for a local maximum distance value nearby the seed point obtained by
the previous method.  
(THIS IS DISABLED AS WITH THE NEW SPHERICITY, THE LOCATIONS OF SPHERICITY MAXIMA
SHOULD MATCH WITH LOCATIONS OF D MAXIMA)
If found within the radius of the distance value at the seed point,
then the seed is relocated to the location of the local maximum.

11/26/08
It returns seeds with highest sphericity as many as 'numseeds.'  When there are multiple
seeds with the exactly same sphericity values, they are all returned.
*/
vector<int>
FindSeedsThresNew(const vector<unsigned char>& L,  //INPUT - foreground segmentation
				  vector<float>& D,          //OUTPUT - distance map
				  vector<float>& Sp,          //OUTPUT - sphericity map
				  vector<unsigned char>& P,  //OUTPUT - local maxima of the sphericity map
				  vector<float>& vSpherical, //OUTPUT - spherical values of seeds
				  int w,
				  int wcore, 
				  float selection_threshold,
				  int numseeds,	
				  const int* dims) {
					  int i;
					  vector<int> vSeedsQ = FindSeedsThres(L, D, Sp, P, vSpherical,
						  w, wcore, selection_threshold, dims);

					  /*if(vSeedsQ.size()>3*numseeds)
					  {
					  vSeedsQ.resize(3*numseeds);
					  }*/
					  float eps = 0.000001;
					  /*i = 0;
					  while(i < numseeds && i<vSpherical.size())
					  {
					  float sval = vSpherical[i];
					  while(i<vSpherical.size() && Abs(vSpherical[i]-sval) < eps)
					  {
					  i++;
					  }
					  //i++;
					  }*/

					  if(vSpherical.empty())
					  {
						  return vSeedsQ;  //failed to pick any seed points.
					  }

					  i = 0; 
					  float thres = Min(0.9, vSpherical[0]);
					  while(i < vSpherical.size() && vSpherical[i] > thres- eps)
					  {
						  i++;
					  }
					  vSeedsQ.resize(3*i);

					  //vSeedsQ = refineSeeds(vSeedsQ, D, L, -1, dims);

					  return vSeedsQ;
}

/*
11/28/08
It uses user supplied click points (default to the center of the sub-volume).

*/
vector<CParticle>
FindSeedsThresNew(const vector<unsigned char>& L,  //INPUT - foreground segmentation
				  vector<float>& D,          //OUTPUT - distance map
				  vector<float>& Sp,          //OUTPUT - sphericity map
				  vector<unsigned char>& P,  //OUTPUT - local maxima of the sphericity map
				  vector<float>& vSpherical, //OUTPUT - spherical values of seeds
				  vector<CParticle>& vclick, //INPUT - click points
				  int w,
				  float selection_threshold,
				  const int* dims) 
{
	int i;
	vector<CParticle> vSeedsQ = FindSeedsThres(L, D, Sp, P, vSpherical, vclick,
		w, selection_threshold, dims);

	float eps = 0.000001;
	if(vSpherical.empty())
	{
		return vSeedsQ;  //failed to pick any seed points.
	}

	//select those above 0.9
	/*i = 0; 
	float thres = Min(0.9, vSpherical[0]);
	while(i < vSpherical.size() && vSpherical[i] > thres- eps)
	{
	i++;
	}
	vSeedsQ.resize(i);*/

	//vSeedsQ = refineSeeds(vSeedsQ, D, L, -1, dims);

	return vSeedsQ;
}

/*
A version on 10/18/2009.
It decouples distance and sphericity maps.
*/
vector<CParticle>
FindSeedsThresV3(const vector<float>& Sp,          //INPUT - sphericity map
				 const vector<unsigned char>& L,          //INPUT - foreground map
				 const vector<float>& D,          //INPUT - distance map
				 const vector<CParticle>& vclick, //INPUT - click points
				 const vector<float>& v,
				 float search_core,
				 float perc,
				 const int* dims) 
{
	vector<CParticle> vseeds;
	for(int c=0; c<vclick.size(); ++c)
	{
		int x0 = vclick[c].m_X;
		int y0 = vclick[c].m_Y;
		int z0 = vclick[c].m_Z;
		float dval = GetData3(D, x0, y0, z0, dims[0], dims[1], dims[2], (float)0);
		float w = 1.5; //Max(3.0, dval*search_core);
		vector<CParticle> vregion = collectConnectedNeighborRegion(L, vclick[c], v, w, dims);
		printf("FindSeedsThresV3: core width = %f, #region voxels = %d\n", w, vregion.size());
		float maxSph = 0;
		for(int i=0; i<vregion.size(); ++i) 
		{ 
			float sval = GetData3(Sp, vregion[i].m_X, vregion[i].m_Y, vregion[i].m_Z, dims[0], dims[1], dims[2], (float)0);
			if(sval > maxSph)
			{
				maxSph = sval;
			}
		}
		float thres = maxSph * perc;
		//printf("FindSeedsThresV3: max sphericity = %f, thres = %f\n", maxSph, thres);
		for(int i=0; i<vregion.size(); ++i) 
		{ 
			float sval = GetData3(Sp, vregion[i].m_X, vregion[i].m_Y, vregion[i].m_Z, dims[0], dims[1], dims[2], (float)0);
			if(sval > thres)
			{
				if(find(vseeds.begin(), vseeds.end(), vregion[i]) == vseeds.end())
				{
					vseeds.push_back(vregion[i]);
					//printf("FindSeedsThresV3: adding (%d, %d, %d)-%f as seeds\n", vregion[i].m_X, vregion[i].m_Y, vregion[i].m_Z, sval);
				}
			}
		}
	}

	return vseeds;
}
