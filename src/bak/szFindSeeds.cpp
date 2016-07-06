
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
          const int* dims) {

  int i, j, k;
  int ndim=3; //3D is assumed
  DistanceTransformFInv(D, L, DistanceTransformMode_Euclid, ndim, dims);
  //printf("Done with Distance Transform.\n");

  int x0=dims[0]/2, y0=dims[1]/2, z0=dims[2]/2;
  float maxSphericity1 = MaximumSphericityValue(1,1,1);
  for(i=0; i<dims[2]; ++i) { 
    for(j=0; j<dims[1]; ++j) { 
      for(k=0; k<dims[0]; ++k) { 
        float val = EvaluateSphericityValue(D,k,j,i,1,1,1,dims)/maxSphericity1;
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
          if(dval>1 && spval>0.0) //d>1 == cannot be at the boarder of the foreground
          //if(!onSurfaceOnSlice(L,i,j,k,dims) && spval>0.0) //d>1 == cannot be at the boarder of the foreground
		  {
            float val=EvaluateSphericityValue(D,k,j,i,sphNbr[0], sphNbr[1], sphNbr[2],dims);
			//printf("(%d, %d, %d): %f, %f, %f\n", k, j, i, dval, spval, val);
           if(val>0)  //valid value should be positive
			{
              float dChecker = Max(Abs(i-z0),Max(Abs(j-y0),Abs(k-x0)));
              if(dChecker>wcore)  //weight the value if it is too far away from the center
			  {
                val*=(1.-(dChecker-wcore)/wcore);
              }
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
			  //printf("Seed at (%d %d %d) %2.2f\n", x+1, y+1, z+1, vvalues[vvalues.size()-i-1]);
		  }
	  }
  }
  return vsub;
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

void
traceEqual(const vector<float>& D,
		   vector<CParticle>& vparticles,
		   int x, int y, int z,
		   float dval,
		   const int* dims)
{
	float dval1 = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	if(dval1 == dval)
	{
		for(int i=0; i<NumNeighbors; ++i)
		{
			int x2 = x+XOffset[i];
			int y2 = y+YOffset[i];
			int z2 = z+ZOffset[i];
			float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
			if(dval2 == dval)
			{
				CParticle pn(x2, y2, z2);
				if(find(vparticles.begin(), vparticles.end(), pn) == vparticles.end())
				{
					vparticles.push_back(pn);
					traceEqual(D, vparticles, x2, y2, z2, dval, dims);
				}
			}
		}
	}
	return;
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
the previous method.  If found within the radius of the distance value at the seed point,
then the seed is relocated to the location of the local maximum.
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

  if(vSeedsQ.size()>3*numseeds)
  {
	  vSeedsQ.resize(3*numseeds);
  }
  vSeedsQ = refineSeeds(vSeedsQ, D, L, -1, dims);

  return vSeedsQ;
}

vector<int>
AscendToBrightest(const vector<unsigned short>& A, 
				  const vector<float>& D, 
				  const vector<int>& vseeds, 
				  const int* dims)
{
	vector<int> refinedSeeds = vseeds;
	for(int i=0; i<vseeds.size(); i+=3)
	{
		int x0 = vseeds[i];
		int y0 = vseeds[i+1];
		int z0 = vseeds[i+2];
		vector<CParticle> vparticles(1, CParticle(x0, y0, z0));
		float dval = GetData3(D, x0, y0, z0, dims[0], dims[1], dims[2], (float)0);
		traceEqual(D, vparticles, x0, y0, z0, dval, dims);
		if(vparticles.size()>1)
		{
			unsigned short valMax = 0;
			int xf = x0;
			int yf = y0;
			int zf = z0;
			for(int j=0; j<vparticles.size(); ++j)
			{
				int x = vparticles[j].m_X;
				int y = vparticles[j].m_Y;
				int z = vparticles[j].m_Z;
				unsigned short val = GetData3(A, x, y, z, dims[0], dims[1], dims[2], (unsigned short)0);
				if(valMax < val)
				{
					xf = x;
					yf = y;
					zf = z;
					valMax = Max(valMax, val);
				}
			}				
			refinedSeeds[i] = xf;
			refinedSeeds[i+1] = yf;
			refinedSeeds[i+2] = zf;
			printf("Refine II: Seed %d is refined from (%d, %d, %d) to (%d, %d, %d)\n",
				i/3+1, vseeds[i]+1, vseeds[i+1]+1, vseeds[i+2]+1, xf+1, yf+1, zf+1);  
		}
	}
	return refinedSeeds;
}
