#include <szCutOutNodule.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDistanceTransform.h>
#include <szConnectedComponent.h>
#include <szMyNeighborOp.h>
#include <szLocalExtrema.h>
#include <szDefaultParam.h>
//#include <szComputeSpecialSphericity.h>
//#include <szDmapSegmentation.h>
//#include <szIsFunctions.h>
//#include <szFitnessScore.h>
#include <szMiscOperations.h>
#include <mex.h>

//const int DefaultLabelVisited = 255;
const float DefaultContinueThreshold = 1.0;
const float DefaultStopThreshold = 0.0;

enum DownHillCode {JustStop, ColorAndContinue, ColorAndStop};

const int OffsetX[6] = {-1, 1, 0, 0, 0, 0};
const int OffsetY[6] = {0, 0, -1, 1, 0, 0};
const int OffsetZ[6] = {0, 0, 0, 0, -1, 1};

//#define printf(a)

void
selectCluster(vector<unsigned char>& L, 
              int x, int y, int z, 
              const int* dims)
{
  int i;
  vector<int> C(L.size());
  int numClusters = ConnectedComponentAnalysisBigger(C, L, NeighborhoodFour, (unsigned char)0, 3, dims);
  if(numClusters > 1)
  {
    int label = GetData3(C, x, y, z, dims[0], dims[1], dims[2], (int) 1);
    for(i=0; i<L.size(); ++i) 
    {
      if(C[i] != label)
        L[i] = 0;
    }
  }
}

/*
To speed up, only provide volume surface to the convex hull process.
*/
void
preConvexHull(vector<unsigned char>& S,   //OUTPUT - must be initialized to zero
              const vector<unsigned char>& C,   //INPUT - segmentation map
              const int* dims) {  
  int i,j,k;
  /*for(i=0; i<C.size(); ++i)
  {
    if(C[i] == DefaultLabelVisited)
      S[i] = 0;
    else
      S[i] = C[i];
  }*/
  
  for(k=0; k<dims[2]; ++k) 
  {
    for(j=0; j<dims[1]; ++j) 
    {
      for(i=0; i<dims[0]; ++i) 
      {
        if(onSurface3(C,i,j,k,dims)) {
          {
            SetData3(S,i,j,k,dims[0],dims[1],dims[2],(unsigned char)1);
          }
        }
      }
    }
  }
}

/*
Segmentation is the intersection of convex-hull and the original foreground segmentation.
This procedure takes care of the taks.
In a rare case, the cluster is broken up into multiple pieces after the above operation.
The function then choose the largest connected component.
*/
void
postConvexHull(vector<unsigned char>& B,    //OUTPUT - must be initialized to zero
               const vector<unsigned char>& C,    //INPUT - convex-hull map
               const vector<unsigned char>& L,   //INPUT - foreground segmentation map
               const int* dims)             //volume size
{
  int i;
  //int cntHollow = 0;
  //int cnt=0;
  for(i=0; i<B.size(); ++i) 
  {
    if(C[i])
	{
		//cnt++;
		if(L[i])
			B[i] = ForegroundColor;
		else
			B[i] = 0;
	}
    else
      B[i] = 0;
  }
  //printf("postConvexHull: %d hollow pixels are found in %d voxels.\n", cntHollow, cnt);
}

/*vector<int>
TraceNew(vector<unsigned char>& B,
		 const vector<unsigned char>& L, 
		 const vector<float>& D, 
		 vector<int> vInd,
		 bool bStrict,
		 const int* dims)
{
  int i,j,k;
  int maxIter = 0;
  for(i=0; i<vInd.size(); i+=3) 
  {
	  int x, y, z;
	  Ind2Sub(x, y, z, vInd[i], dims);
	  float dVal = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float) 0);
	  maxIter = Max(maxIter, int(dVal-0.0001)+1);
	  SetData3(B, x, y, z, dims[0], dims[1], dims[2], (unsigned char)1);
	  //printf("TraceNew: (%d %d %d) %f\n", x+1, y+1, z+1, dVal);
  }  
  //printf("TraceNew: Max iteration = %d (%d)\n", maxIter, vInd.size());
  vector<int> vNInd = MakeFourNeighborhood(3, dims);

  int iter = 0;
  while(vInd.size()>0 && iter++ <= maxIter)
  {
    vector<int> vInd2 = vInd;
    vInd.clear();
    for(i=0; i<vInd2.size(); ++i)
    {
      int x2, y2, z2;
      Ind2Sub(x2, y2, z2, vInd2[i], dims);
      float dVal1 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float) 0);

      for(j=0; j<vNInd.size(); ++j)
      {
        int ind = vInd2[i] + vNInd[j];
        if(ind>= 0 && ind<B.size())
        {
          if(B[ind] == 0 && L[ind]) 
          {
            int x3, y3, z3;
            Ind2Sub(x3, y3, z3, ind, dims);
            float dVal2 = GetData3(D, x3, y3, z3, dims[0], dims[1], dims[2], (float) 0);
            //printf("(%d %d %d: %2.2f) -> (%d %d %d: %2.2f):\t",
            //  x2+1, y2+1, z2+1, dVal1, x3+1, y3+1, z3+1, dVal2);
            //provide two different threshold for dVal>1 and dVal==1, due to 
            //inaccurate gradient at the boundary of the foreground.
            if((bStrict && dVal2 < dVal1) || (!bStrict && dVal2 <= dVal1))
            {
				//printf("IN:\n");
				//unsigned char bVal = GetData3(B, x3, y3, z3, dims[0], dims[1], dims[2], (unsigned char) 0);
				SetData3(B, x3, y3, z3, dims[0], dims[1], dims[2], (unsigned char)1);
				if(find(vInd.begin(), vInd.end(), ind)==vInd.end())
				{
					vInd.push_back(ind);
				}
            }
            else
            {
              //printf("OUT:\n");
            }
          }
        }
      }
    }
  }
  return vInd;
}*/

bool
AlwaysTrace(float a, float b)
{
	return true;
}

bool
StrictDownHillTrace(float a, float b)
{
	return a > b ? true: false;
}

bool
NonStrictDownHillTrace(float a, float b)
{
	return a >= b ? true: false;
}

/*vector<CParticle>
ColorBall(vector<unsigned char>& B,
		  const vector<unsigned char>& L, 
		  int x, int y, int z,
		  float radius,
		  int ndim,
		  const int* dims)
{
	vector<unsigned char> C(B.size());
	vector<CParticle> vseeds;
	float eps = 0.1;
	int i;
	if(SetData3(C, x, y, z, dims[0], dims[1], dims[2], (unsigned char)1))
	{
		vector<float> M(B.size());
		DistanceTransformF(M, C, DistanceTransformMode_Euclid, ndim, dims);
		for(i=0; i<B.size(); ++i)
		{
			if(M[i] <= radius && L[i])
			{
				B[i] = 1;
				if(M[i] < radius+eps && M[i] > radius-eps)
				{
					int x2, y2, z2;
					Ind2Sub(x2, y2, z2, i, dims);
					CParticle sd(x2, y2, z2, 0);
					vseeds.push_back(sd);
				}
			}
		}
	}
	return vseeds;
}*/

/*
This version set a life time for each trace.  The life time is dependent on the distance value at the
seed of the trace.
*/
/*vector<CParticle>
TraceNew(vector<unsigned char>& B,
		 const vector<unsigned char>& L, 
		 const vector<float>& D, 
		 vector<CParticle> vseeds,
		 float scale,
		 bool (*decision_func)(float, float),
		 const int* dims)
{
	int i,j,k;
	int maxIter = 0;
	for(i=0; i<vseeds.size(); i++) 
	{
		int x = vseeds[i].m_X;
		int y = vseeds[i].m_Y;
		int z = vseeds[i].m_Z;
		float dVal = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float) 0);
		if(dVal > 0)
		{
			vseeds[i].m_Life = (int)(scale * dVal);
			SetData3(B, x, y, z, dims[0], dims[1], dims[2], (unsigned char)1);
		}
		else
		{
			vseeds[i].m_Life = 0;
		}
	}  
	//const int XOffset[] = {-1, 0, 1, -1, 1, -1, 0, 1, 0, 0};
	//const int YOffset[] = {-1, -1, -1, 0, 0, 1, 1, 1, 0, 0};
	//const int ZOffset[] = {0, 0, 0, 0, 0, 0, 0, 0, -1, 1};
	const int XOffset[] = {-1, 0, 1, 0, 0, 0};
	const int YOffset[] = {0, -1, 0, 1, 0, 0};
	const int ZOffset[] = {0, 0, 0, 0, -1, 1};
	const int NumNeighbors = sizeof(XOffset)/sizeof(XOffset[0]);

	vector<CParticle> vdied;
	while(vseeds.size()>0)
	{
		vector<CParticle> vseeds2 = vseeds;
		vseeds.clear();
		for(i=0; i<vseeds2.size(); ++i)
		{
			int x2 = vseeds2[i].m_X;
			int y2 = vseeds2[i].m_Y;
			int z2 = vseeds2[i].m_Z;
			float dVal1 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float) 0);
			for(j=0; j<NumNeighbors; ++j)
			{
				int x3 = x2 + XOffset[j];
				int y3 = y2 + YOffset[j];
				int z3 = z2 + ZOffset[j];
				unsigned char lVal2 = GetData3(L, x3, y3, z3, dims[0], dims[1], dims[2], (unsigned char)0);
				if(lVal2) 
				{
					float dVal2 = GetData3(D, x3, y3, z3, dims[0], dims[1], dims[2], (float) 0);
					if(decision_func(dVal1, dVal2))
					{
						SetData3(B, x3, y3, z3, dims[0], dims[1], dims[2], (unsigned char)ForegroundColor);
						if(vseeds2[i].m_Life > 1)
						{
							CParticle sd(x3, y3, z3, vseeds2[i].m_Life-1);
							bool found = false;
							for(k=0; k<vseeds.size(); ++k)
							{
								if(sd == vseeds[k])
								{
									vseeds[k].m_Life = Max(sd.m_Life, vseeds[k].m_Life);
									found = true;
								}
							}
							if(!found)
							{
								vseeds.push_back(sd);
							}
						}
						else
						{
							CParticle sd(x3, y3, z3, vseeds2[i].m_Life-1);
							bool found = false;
							for(k=0; k<vdied.size(); ++k)
							{
								if(sd == vdied[k])
								{
									vdied[k].m_Life = Max(sd.m_Life, vdied[k].m_Life);
									found = true;
								}
							}
							if(!found)
							{
								vdied.push_back(sd);
							}
						}
					}
					else
					{
					}
				}
			}
		}
	}
	return vdied;
}*/

/*
A slight modification to TraceNew.  Here, only those traces whose distance values are less than 
certain threshold derived from the data.  
We set the threshold to the 1st quartile of the distance values at seeds.  The motivation is to
prevent over-segmentation of wall-attached nodules by only allowing those seeds that can trace 
radially towards the boundary.
The rationale of the threshold is that if a 1/4 of the boundary is properly segmented, then the
segmentation will be considered correct.
*/
vector<CParticle>
TraceNew2(vector<unsigned char>& B,
		 const vector<unsigned char>& L, 
		 const vector<float>& D, 
		 vector<CParticle> vseeds,
		 float scale,
		 bool (*decision_func)(float, float),
		 const int* dims)
{
	int i,j,k;
	int maxIter = 0;
	vector<float> vDvalues;
	for(i=0; i<vseeds.size(); i++) 
	{
		int x = vseeds[i].m_X;
		int y = vseeds[i].m_Y;
		int z = vseeds[i].m_Z;
		float dVal = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float) 0);
		vDvalues.push_back(dVal);
	}
	sort(vDvalues.begin(), vDvalues.end());
	float threshold = vDvalues[vDvalues.size()/4]; //use 1st quartile to prevent wall-leak
	int cnt = 0;
	for(i=0; i<vseeds.size(); i++) 
	{
		int x = vseeds[i].m_X;
		int y = vseeds[i].m_Y;
		int z = vseeds[i].m_Z;
		float dVal = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float) 0);
		if(dVal >= threshold)
		{
			//vseeds[i].m_Life = Max((int)ceil(scale * dVal), 2);
			vseeds[i].m_Life = Round(scale * dVal);
			vseeds[i].m_Life = Max(1.01, vseeds[i].m_Life);
			vseeds[i].m_X0 = x;
			vseeds[i].m_Y0 = y;
			vseeds[i].m_Z0 = z;
			SetData3(B, x, y, z, dims[0], dims[1], dims[2], (unsigned char)1);
			cnt++;
		}
		else
		{
			vseeds[i].m_Life = 0;
		}
	}  
	//printf("TraceNew2: threshold = %f, #seeds = %d\n", threshold, cnt);
	const int XOffset[] = {-1, 0, 1, 0, 0, 0};
	const int YOffset[] = {0, -1, 0, 1, 0, 0};
	const int ZOffset[] = {0, 0, 0, 0, -1, 1};
	const int NumNeighbors = sizeof(XOffset)/sizeof(XOffset[0]);

	vector<CParticle> vdied;
	while(vseeds.size()>0)
	{
		vector<CParticle> vseeds2 = vseeds;
		vseeds.clear();
		for(i=0; i<vseeds2.size(); ++i)
		{
			int x2 = vseeds2[i].m_X;
			int y2 = vseeds2[i].m_Y;
			int z2 = vseeds2[i].m_Z;
			int x0 = vseeds2[i].m_X0;
			int y0 = vseeds2[i].m_Y0;
			int z0 = vseeds2[i].m_Z0;
			float dVal1 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float) 0);
			bool bDied = true;
			for(j=0; j<NumNeighbors; ++j)
			{
				int x3 = x2 + XOffset[j];
				int y3 = y2 + YOffset[j];
				int z3 = z2 + ZOffset[j];
				if((x3-x0)*(x3-x0)+(y3-y0)*(y3-y0)+(z3-z0)*(z3-z0)>=vseeds2[i].m_Life*vseeds2[i].m_Life)
				{
					continue;
				}
				unsigned char lVal2 = GetData3(L, x3, y3, z3, dims[0], dims[1], dims[2], (unsigned char)0);
				unsigned char bVal2 = GetData3(B, x3, y3, z3, dims[0], dims[1], dims[2], (unsigned char)0);
				if(lVal2 && !bVal2) 
				{
					float dVal2 = GetData3(D, x3, y3, z3, dims[0], dims[1], dims[2], (float) 0);
					if(decision_func(dVal1, dVal2))
					{
						SetData3(B, x3, y3, z3, dims[0], dims[1], dims[2], (unsigned char)ForegroundColor);
						CParticle sd(x3, y3, z3, vseeds2[i].m_Life, x0, y0, z0);
						bool found = false;
						for(k=0; k<vseeds.size(); ++k)
						{
							if(sd == vseeds[k])
							{
								found = true;
								if(sd.m_Life > vseeds[k].m_Life)
								{
									vseeds[k] = sd;
								}
							}
						}
						if(!found)
						{
							vseeds.push_back(sd);
						}
						bDied = false;
					}
				}
			}
			if(bDied)
			{
				CParticle sd(x2, y2, z2, (int)(scale * dVal1), x2, y2, z2);
				vdied.push_back(sd);
			}
		}
	}
	return vdied;
}

void
MergeLabels(vector<unsigned char>& A,
            vector<unsigned char>& B,
            const int* dims)
{
  for(int i=0; i<A.size(); ++i)
  {
    if(B[i])
    {
      A[i] = ForegroundColor;
    }
  }
}


