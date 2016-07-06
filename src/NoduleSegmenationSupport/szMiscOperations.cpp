#ifdef MEX_DLL
#include <mex.h>
#endif
#include <math.h>
#include <algorithm>
using namespace std;

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szConnectedComponent.h>
#include <szConnectedComponentC.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szEigen.h>

int GetOffsetX6(int n) 
{
  switch(n) {
  case 1:
    return -1;
  case 2:
    return 1;
  case 0:
  case 3:
  case 4:
  case 5:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetY6(int n) 
{
  switch(n) {
  case 0:
    return -1;
  case 3:
    return 1;
  case 1:
  case 2:
  case 4:
  case 5:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetZ6(int n) 
{
  switch(n) {
  case 4:
    return -1;
  case 5:
    return 1;
  case 0:
  case 1:
  case 2:
  case 3:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
  //const int zoff[6]={0,0,0,0,-1,1};
  //return zoff[n];
}

int GetOffsetX10(int n) 
{
  switch(n) {
  case 1:
  case 4:
  case 6:
    return -1;
  case 3:
  case 5:
  case 8:
    return 1;
  case 0:
  case 2:
  case 7:
  case 9:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetY10(int n) 
{
  switch(n) {
  case 1:
  case 2:
  case 3:
    return -1;
  case 6:
  case 7:
  case 8:
    return 1;
  case 0:
  case 4:
  case 5:
  case 9:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
}
 
int GetOffsetZ10(int n) 
{
  switch(n) {
  case 0:
    return -1;
  case 9:
    return 1;
  case 1:
  case 2:
  case 3:
  case 4:
  case 5:
  case 6:
  case 7:
  case 8:
    return 0;
  default:
    //mexPrintf("Illegal neighbor index of %d for 6-neighborhood\n");
    return 0;
  }
  //const int zoff[6]={0,0,0,0,-1,1};
  //return zoff[n];
}

int GetOffsetX26(int n) 
{
  if(n<9)
    return (n%3)-1;
  else if(n<17) {
    switch(n) {
    case 9:
    case 12:
    case 14:
      return -1;
    case 10:
    case 15:
      return 0;
    case 11:
    case 13:
    case 16:
      return 1;
	default:
		return 0;
    }
  }
  else 
    return (n-17)%3 - 1;
}
 
int GetOffsetY26(int n) 
{
  if(n<9)
    return (n/3)-1;
  else if(n<17) {
    switch(n) {
    case 9:
    case 10:
    case 11:
      return -1;
    case 12:
    case 13:
      return 0;
    case 14:
    case 15:
    case 16:
      return 1;
	default:
		return 0;
    }
  }
  else 
    return (n-17)/3-1;
}
 
int GetOffsetZ26(int n) {
  if(n<9)
    return -1;
  else if(n<17)
    return 0;
  else //if(n<26)
    return 1;
}


/*
sort values and returns sorted indices.
*/
vector<int>
indexedSort(vector<float>& vvalues) 
{
  vector<indexedData> vdata(vvalues.size());
  int i;
  for(i=0; i<vvalues.size(); ++i) {
    vdata[i].data = vvalues[i];
    vdata[i].index = i;
  }

  sort(vdata.begin(),vdata.end());

  vector<int> vindex(vvalues.size());
  for(i=0; i<vvalues.size(); ++i) {
    vvalues[i] = vdata[i].data;
    vindex[i] = vdata[i].index;
  }

  return vindex;
}

/*
Check if the voxel is on the boundary of the segmentation.
*/

bool
onSurface3(const vector<unsigned char>& A, int x, int y, int z, const int* dims) 
{
  if(GetData3(A,x,y,z,dims[0],dims[1],dims[2],(unsigned char)0)) {
    if(!GetData3(A,x-1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(!GetData3(A,x+1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(!GetData3(A,x,y-1,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(!GetData3(A,x,y+1,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(!GetData3(A,x,y,z-1,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(!GetData3(A,x,y,z+1,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
  }
  return false;
}

bool
onSurface3(const vector<int>& A, int x, int y, int z, const int* dims) 
{
  if(GetData3(A,x,y,z,dims[0],dims[1],dims[2],(int)0)) {
    if(!GetData3(A,x-1,y,z,dims[0],dims[1],dims[2],(int)0))
      return true;
    if(!GetData3(A,x+1,y,z,dims[0],dims[1],dims[2],(int)0))
      return true;
    if(!GetData3(A,x,y-1,z,dims[0],dims[1],dims[2],(int)0))
      return true;
    if(!GetData3(A,x,y+1,z,dims[0],dims[1],dims[2],(int)0))
      return true;
    if(!GetData3(A,x,y,z-1,dims[0],dims[1],dims[2],(int)0))
      return true;
    if(!GetData3(A,x,y,z+1,dims[0],dims[1],dims[2],(int)0))
      return true;
  }
  return false;
}

bool
onBoundary3(const vector<unsigned char>& A, int x, int y, int z, const int* dims) 
{
  if(!GetData3(A,x,y,z,dims[0],dims[1],dims[2],(unsigned char)1)) {
    if(GetData3(A,x-1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x+1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y-1,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y+1,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y,z-1,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y,z+1,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
  }
  return false;
}

/*
Check if the voxel is on the boundary of the segmentation.
*/
bool
onSurface3(const vector<unsigned char>& A, int id, const int* dims) 
{
  int ndim=3;
  vector<int> vsub = Ind2Sub(id,ndim,dims);
  if(!GetData3(A,vsub[0]-1,vsub[1],vsub[2],dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,vsub[0]+1,vsub[1],vsub[2],dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,vsub[0],vsub[1]-1,vsub[2],dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,vsub[0],vsub[1]+1,vsub[2],dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,vsub[0],vsub[1],vsub[2]-1,dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,vsub[0],vsub[1],vsub[2]+1,dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  
  return false;
}

/*
Check if the background voxel is touching the foreground.
*/
bool
onAdjacentSurface3(const vector<unsigned char>& A, int x, int y, int z, const int* dims) 
{
  if(GetData3(A,x,y,z,dims[0],dims[1],dims[2],(unsigned char)1) == 0) {
    if(GetData3(A,x-1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x+1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y-1,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y+1,z,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y,z-1,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
    if(GetData3(A,x,y,z+1,dims[0],dims[1],dims[2],(unsigned char)0))
      return true;
  }
  return false;
}

/*
Check if the voxel is on the boundary of the segmentation within the current
axial slice.
*/
bool
onSurfaceOnSlice(const vector<unsigned char>& A, int x, int y, int z, const int* dims) {
  if(!GetData3(A,x,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
    return false;
  if(!GetData3(A,x-1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,x+1,y,z,dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,x,y-1,z,dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  if(!GetData3(A,x,y+1,z,dims[0],dims[1],dims[2],(unsigned char)0))
    return true;
  
  return false;
}

/*
Compute the diameter of the segmentation in the current axial slice.
This version uses an ELCAP protocol to approximate the radii.
*/
int
ComputeAxialDiameterDirect(double& major,            //OUTPUT - major diameter
                           double& minor,            //OUTPUT - minor diameter
                           const vector<int>& vx,           //x coordinates
                           const vector<int>& vy,           //y cordinates
                           const vector<float>& v) 
{
  //printf("Enter ComputeAxialDiameter()\n");

  if(vx.empty()) {
    //printf("voxel coordinate empty.  Cannot continue.\n");
    major = 0;
    minor = 0;
    return 0;
  }
  else if(vx.size()==1) {
    //printf("only single point is in the cluster.\n");
    major = 1.0;
    minor = 1.0;
    return 1;
  }
  else if(vx.size()!=vy.size()) {
    //printf("voxel coordinate array mismatch.  Cannot continue.\n");
    major = 1.0;
    minor = 1.0;
    return -1;
  }
  else {
    //first compute the major axis - brute force
    //compare every pair and get the one that is furtherest apart
    double max_dist = 0;
    int id1, id2;
    double dX, dY;
    int X1, Y1;
    bool found = false;
	int i;
    for(i=0; i<vx.size()-1; ++i) 
    {
      int x1=vx[i];
      int y1=vy[i]; 
      for(int j=i+1; j<vx.size(); ++j) 
      {
        double dx = v[0]*(x1-vx[j]);
        double dy = v[1]*(y1-vy[j]);
        double dist = (double) dx*dx+dy*dy;
        if(dist>max_dist) {
          id1 = i;
          id2 = j;
          dX = dx;
          dY = dy;
          X1 = x1;
          Y1 = y1;
          max_dist = dist;
          found = true;
        }
      }
    }
    //printf("Maximum distance computed\n");
    if(!found || max_dist<=0) { //this should never happen
      //printf("No axial diameter found...\n");
      major=1.0;
      minor=1.0;
      return -1;
    }

    major = sqrt(max_dist);
    dX /= major;
    dY /= major;
    //printf("Major axis computed\n");

    //next compute the minor axis
    //use the major axix, find a point on each side of the axis, 
    //whose projection length (i.e. perp to the major axis) is longest
    double max_dist2 = 0;
    double min_dist2 = 0;
    bool found1 = false;
    bool found2 = false;
    for(i=0; i<vx.size(); ++i) 
    {
      if(i==id1 || i==id2)
        continue;

      double dx = v[0]*(vx[i]-X1);
      double dy = v[1]*(vy[i]-Y1);
      double dist = -dx*dY + dy*dX; //the length of projection to the major segment
      if(dist>max_dist2) 
        max_dist2 = dist;
      else if(dist<min_dist2)
        min_dist2 = dist;
    }

	//add 1 so that a small nodule of width 1 results in minor length of 1.
    minor = max_dist2 - min_dist2 + 1;
	major += 1.0; 
	//printf("ComputeAxialDiameterDirect: major = %f, minor = %f\n", major, minor);
    
    //printf("Exit ComputeAxialDiameter()\n");
    return 1;
  }
}

/*
Compute the diameter of the segmentation in the current axial slice.
This version uses Gaussian ellipsoid fit to approximate the radii.
*/
int
ComputeAxialDiameterEllipse(double& major,            //OUTPUT - major diameter
                            double& minor,            //OUTPUT - minor diameter
                            const vector<int>& vx,           //x coordinates
                            const vector<int>& vy,           //y cordinates
                            const vector<float>& v) 
{
  //printf("Enter ComputeAxialDiameter()\n");

  if(vx.empty()) 
  {
    //printf("voxel coordinate empty.  Cannot continue.\n");
    major = 0;
    minor = 0;
    return 0;
  }
  else if(vx.size()==1) 
  {
    //printf("only single point is in the cluster.\n");
    major = 1.0;
    minor = 1.0;
    return 1;
  }
  else if(vx.size()!=vy.size()) 
  {
    //printf("voxel coordinate array mismatch.  Cannot continue.\n");
    major = 1.0;
    minor = 1.0;
    return -1;
  }
  else 
  {
    double sx=0, sy=0, sxx=0, syy=0, sxy=0;
    for(int i=0; i<vx.size(); ++i) 
    {
      double x=v[0] * vx[i];
      double y=v[1] * vy[i];
      sx+=x;
      sxx+=x*x;
      sy+=y;
      syy+=y*y;
      sxy+=x*y;
    }
    sx/=vx.size();
    sy/=vy.size();
    sxx=sxx/vx.size()-sx*sx;
    syy=syy/vy.size()-sy*sy;
    sxy=sxy/vx.size()-sx*sy;

    //double ee = sqrt((sxx-syy)*(sxx-syy)+sxy*sxy);
    //major = 3.2*sqrt((sxx+syy)+ee);
    //minor = 3.2*sqrt((sxx+syy)-ee);
    double ee = sqrt((sxx-syy)*(sxx-syy)+4.0*sxy*sxy);
    major = sqrt(0.5*((sxx+syy)+ee));
    minor = sqrt(0.5*((sxx+syy)-ee));
    return 1;
  }
}

/*
*/
float
EvaluateSphericityValue(const vector<float>& D,
                        int x,
                        int y, 
                        int z,
                        int widthX,
                        int widthY,
                        int widthZ,
                        const int* dims) 
{

  bool b;
  float v0 = GetData3(D,x,y,z,dims[0],dims[1],dims[2],b);
  if(b==false && v0<1) //disregard the background pixel
    return 0;

  float sum=0;
  int cntFG=0; //foreground neighbor counter
  for(int i=z-widthZ; i<=z+widthZ; ++i) 
  {
    if(i>=0 && i<dims[2]) {
      for(int j=y-widthY; j<=y+widthY; ++j) 
      {
        if(j>=0 && j<dims[1]) {
          for(int k=x-widthX; k<=x+widthX; ++k) 
          {
            if(k>=0 && k<dims[0]) 
            {
              float v2 = GetData3(D,k,j,i,dims[0],dims[1],dims[2],b);
              if(b && v2>0) //added v2>0 condition on 01/06/09 
              {
				  if(v0 > v2)
					  sum += 1.0;
				 // else if(v0 < v2)
					//  sum -= 1.0;
              }
            }
          }
        }
      }
    }
  }

  return sum;
}

/*
Dot product based measure.
*/
float
EvaluateSphericityValueOrig(const vector<float>& D,
                        int x,
                        int y, 
                        int z,
                        int widthX,
                        int widthY,
                        int widthZ,
                        const int* dims) 
{

  bool b;
  float v0 = GetData3(D,x,y,z,dims[0],dims[1],dims[2],b);
  if(b==false && v0<1) //disregard the background pixel
    return 0;

  float sum=0;
  int cntFG=0; //foreground neighbor counter
  for(int i=z-widthZ; i<=z+widthZ; ++i) 
  {
    if(i>=0 && i<dims[2]) {
      for(int j=y-widthY; j<=y+widthY; ++j) 
      {
        if(j>=0 && j<dims[1]) {
          for(int k=x-widthX; k<=x+widthX; ++k) 
          {
            if(k>=0 && k<dims[0]) 
            {
              float v2 = GetData3(D,k,j,i,dims[0],dims[1],dims[2],b);
              if(b) 
              {
                //if(v2>0) {
                  float dd=(i-z)*(i-z)+(j-y)*(j-y)+(k-x)*(k-x);
                  float df=v0-v2;
                  sum += df*sqrt(dd);
                //}
                /*else if(v2==0 && (Abs(i-z)==1 || Abs(j-y)==1 || Abs(k=x)==1)) {
                  float dd=(i-z)*(i-z)+(j-y)*(j-y)+(k-x)*(k-x);
                  float df=v0-v2;
                  sum += df*sqrt(dd);
                }*/
              }
            }
          }
        }
      }
    }
  }

  return sum;
}

/*
*/
float
MaximumSphericityValue(int widthX,
                       int widthY,
                       int widthZ) 
{
	return (2*widthX+1)*(2*widthY+1)*(2*widthZ+1) - 1;
}

/*
Dot product based measure.
*/
float
MaximumSphericityValueOrig(int widthX,
                       int widthY,
                       int widthZ) 
{
  
  float sum=0;
  for(int i=-widthZ; i<=widthZ; ++i) 
  {
    for(int j=-widthY; j<=+widthY; ++j) 
    {
      for(int k=-widthX; k<=+widthX; ++k) 
      {
        {
          float dd=i*i+j*j+k*k;
          sum += dd;
        }
      }
    }
  }
  return sum;
}

/*
select a connected component that includes (x, y, z).
*/
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
From a set of positions, it returns a representative position defined as:
the position in the set that is closest to the centroid of the positions in
the set.
*/

vector<int>
representativePosition(const vector<int>& pos) 
{
  vector<int> vrep;

  if(pos.size()<3) 
  {
    //needs at least one position (x,y,z)
    //return empty vector
    return vrep;
  }
  else if(pos.size()<6) 
  {
    //single position - this is the representative one
    return pos;
  }

  double ax=0;
  double ay=0;
  double az=0;
  int count=0;
  int n;
  for(n=0; n<pos.size(); n+=3) 
  {
    ax+=pos[n];
    ay+=pos[n+1];
    az+=pos[n+2];
    count++;
  }
  ax/=count;
  ay/=count;
  az/=count;

  int rx=pos[0];
  int ry=pos[1];
  int rz=pos[2];
  double dst=(ax-rx)*(ax-rx)+(ay-ry)*(ay-ry)+(az-rz)*(az-rz);
  for(n=3; n<pos.size(); n+=3) 
  {
    int x=pos[n];
    int y=pos[n+1];
    int z=pos[n+2];
    double dst2=(ax-x)*(ax-x)+(ay-y)*(ay-y)+(az-z)*(az-z);
    if(dst2<dst) 
    {
      rx=x;
      ry=y;
      rz=z;
    }
  }

  vrep.push_back(rx);  
  vrep.push_back(ry);  
  vrep.push_back(rz);

  return vrep;
}


/*
The routine tries to fill in missed part of the foreground near the boundary due to 
inmature termination of the region growing.
*/
vector<unsigned char>
refineSegmentation(vector<unsigned char>& S, //segmentation
				   const vector<unsigned char>& L, //foreground
				   const vector<int>& cn, //a causal neighborhood as a structural element
				   int ndim, 
				   const int* dims //dimension of the volume
				   )
{
	int nvoxels = numberOfElements(ndim,dims);
	int i, j, k;
	vector<unsigned char> B = S;
	//printf("entering refineSegmentation.\n");
	int cnt = 0;
	for(i=0; i<nvoxels; ++i)
	{
		if(!S[i] && L[i])
		{
			for(int m=0; m<cn.size(); ++m)
			{
				if(NeighborCheck(i, cn[m], ndim, dims) && NeighborCheck(i, -cn[m], ndim, dims))
				{
					if((!L[i+cn[m]] && S[i-cn[m]]) || (!L[i-cn[m]] && S[i+cn[m]]))
					{
						B[i] = ForegroundColor;
						cnt++;
						//printf("refineSegmentation: paint %d\n", i);
					}
				}
			}
		}
	}

	return B;
}

/*float
computeForegroundDScore(const vector<int>& vseeds,
						const int* dims)
{
	float score = 0;
	int cnt = 0;
	for(int i=0; i<vseeds.size(); i+=3)
	{
		int x = vseeds[i];
		int y = vseeds[i+1];
		int z = vseeds[i+2];
		float dd = (x-dims[0]*.5)*(x-dims[0]*.5)+(y-dims[1]*.5)*(y-dims[1]*.5)+(z-dims[2]*.5)*(z-dims[2]*.5);
		score += sqrt(dd);
		cnt ++;
	}
	if(cnt)
	{
		score /= (float)cnt;
		float wgt = 0.25;
		return 1.0/(wgt*score+1);
	}
	else
		return 0;
}*/

float
computeForegroundSScore(const vector<float>& vSp)
{
	float score;
	if(vSp.empty())
		score = 0;
	else
		score = vSp[0];
	return score;
}

/*int 
compareForegroundFitness(const vector<float>& vSp1,
						 const vector<int>& vseeds1,
						 const vector<float>& vSp2,
						 const vector<int>& vseeds2,
						 const int* dims)
{
	float dscore1 = computeForegroundDScore(vseeds1, dims);
	float sscore1 = computeForegroundSScore(vSp1);
	float dscore2 = computeForegroundDScore(vseeds2, dims);
	float sscore2 = computeForegroundSScore(vSp2);
	printf("compareForegroundFitness: %f vs %f (%f+%f vs %f+%f)(D+S)\n", dscore1+sscore1,dscore2+sscore2, dscore1,sscore1,dscore2,sscore2);
	if(dscore1 + sscore1 > dscore2 + sscore2)
	{
		return 1;
	}
	else if(dscore1 + sscore1 < dscore2 + sscore2)
	{
		return -1;
	}
	else
	{
		if(dscore2 + sscore2 > 0)
		{
			return -1;
		}
		else
		{
			return 0;
		}
	}
}*/

float
boostForegroundScore(float score, int noduleType)
{
	float boostSolidCoeff = 1.2;
	float boostNonSolidCoeff = 1.1;
	switch(noduleType)
	{
	case SZ_LABEL_SOLID_NODULE:
			/*if(score>0.6)
			{
				score = 2.0;
			}
			else
			{
				score *= boostSolidCoeff;
			}*/
		score *= boostSolidCoeff;
		break;
	case SZ_LABEL_NONSOLID_NODULE:
		/*if(score>0.6)
		{
		score = 1.0;
		}
		else
		{
		score *= boostNonSolidCoeff;
		}
		break;*/
		score *= boostNonSolidCoeff;
		break;
	}
	return score;
};

/*
compare the sphericity values.
incorporate priority of solid over non-solid and part-solid,
and priority of non-solid over part-solid
*/
int 
compareForegroundFitness2(const vector<float>& vSp1,
						  int noduleType1,
						  const vector<float>& vSp2,
						  int noduleType2)
{
	float sscore1 = computeForegroundSScore(vSp1);
	float sscoreB1 = boostForegroundScore(sscore1, noduleType1);

	float sscore2 = computeForegroundSScore(vSp2);
	float sscoreB2 = boostForegroundScore(sscore2, noduleType2);

	printf("compareForegroundFitness2: %f (%f) vs %f (%f)\n", sscoreB1, sscore1, sscoreB2, sscore2);
	if(sscoreB1 > sscoreB2)
	{
		return 1;
	}
	else if(sscoreB1 < sscoreB2)
	{
		return -1;
	}
	else
	{
		if(sscoreB1 > 0) //tie but non-zero
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
}

vector<double>
computeCovariance(const vector<double>& vx,
				  const vector<double>& vy,
				  const vector<double>& vz)
{
	int i;
	double mx=0, my=0, mz=0;
	double exx=0, eyy=0, ezz=0, exy=0, exz=0, eyz=0;
	for(i=0; i<vx.size(); ++i)
	{
		mx += vx[i];
		my += vy[i];
		mz += vz[i];
		exx += vx[i]*vx[i];
		exy += vx[i]*vy[i];
		exz += vx[i]*vz[i];
		eyy += vy[i]*vy[i];
		eyz += vy[i]*vz[i];
		ezz += vz[i]*vz[i];
	}
	mx = mx/vx.size();
	my = my/vx.size();
	mz = mz/vx.size();
	exx = exx/vx.size() - mx*mx;
	exy = exy/vx.size() - mx*my;
	exz = exz/vx.size() - mx*mz;
	eyy = eyy/vx.size() - my*my;
	eyz = eyz/vx.size() - my*mz;
	ezz = ezz/vx.size() - mz*mz;

	vector<double> vCov(9);
	vCov[0] = exx; 
	vCov[1] = exy; 
	vCov[2] = exz; 
	vCov[3] = exy; 
	vCov[4] = eyy; 
	vCov[5] = eyz; 
	vCov[6] = exz; 
	vCov[7] = eyz; 
	vCov[8] = ezz; 

	return vCov;
}

/*
Compute the center of mass of a segmentation volume.
*/
vector<double>
computeCentroid(const vector<unsigned char>& S, //segmentation
				int ndim, //should be 3
				const int* dims)
{
	vector<double> centroid(3,0);
	int i, j, k;
	int cnt = 0;
	for(k=0; k<dims[2]; ++k) {
		for(i=0; i<dims[1]; ++i) {
			for(j=0; j<dims[0]; ++j) {
				if(GetData3(S, j, i, k, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					centroid[0] += j;
					centroid[1] += i;
					centroid[2] += k;
					cnt++;
				}
			}
		}
	}

	if(cnt)
	{
		centroid[0] /= (double) cnt;
		centroid[1] /= (double) cnt;
		centroid[2] /= (double) cnt;
	}
	return centroid;
}

/*
a small hole can alter the distance transform.  So plung them.
*/
void
FillBackgroundHoles(vector<unsigned char>& B,
					int minSize,
					int label,
					const int* dims)
{
	int imageSize = dims[0]*dims[1];
	vector<unsigned char> BP(imageSize, 0);
	vector<int> CP(imageSize, 0);
	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<imageSize; ++j)
		{
			BP[j] = B[i*imageSize + j];
			CP[j] = 0;
		}
		vector<int> vcount;
		int numClusters = ConnectedComponentAnalysisEqual(CP, BP, vcount, NeighborhoodFour, (unsigned char)0, 2, dims);
		//int numClusters = ConnectedComponentAnalysisEqual(CP, BP, NeighborhoodFour, (unsigned char)0, 2, dims);

		vector<int> bFill(numClusters+1, true); //clusters that touch image boundaries will not be filled.
		//check the four image boundary and mark with ones should be ignored.
		for(int k=0; k<dims[1]; ++k)
		{
			int m = GetData2(CP, 0, k, dims[0], dims[1], 0);
			if(m > 0) bFill[m] = false;
			m = GetData2(CP, dims[0]-1, k, dims[0], dims[1], 0);
			if(m > 0) bFill[m] = false;
		}
		for(int k=0; k<dims[0]; ++k)
		{
			int m = GetData2(CP, k, 0, dims[0], dims[1], 0);
			if(m > 0) bFill[m] = false;
			m = GetData2(CP, k, dims[1]-1, dims[0], dims[1], 0);
			if(m > 0) bFill[m] = false;
		}
		for(int j=0; j<imageSize; ++j)
		{
			if(CP[j] && bFill[CP[j]] && vcount[CP[j]] < minSize)
			{
				B[i*imageSize + j] = label;
			}
		}
	}
}

float
trilinearInterpolation(const vector<float>& D,
					   float x, 
					   float y, 
					   float z,
					   const int* dims)
{
	int x0 = (int)x;
	int y0 = (int)y;
	int z0 = (int)z;
	int x1 = x0+1;
	int y1 = y0+1;
	int z1 = z0+1;
	float r = x-x0;
	float s = y-y0;
	float t = z-z0;

	float dvals[8];
	dvals[0] = GetData3(D, x0, y0, z0, dims[0], dims[1], dims[2], (float)0);
	dvals[1] = GetData3(D, x1, y0, z0, dims[0], dims[1], dims[2], (float)0);
	dvals[2] = GetData3(D, x0, y1, z0, dims[0], dims[1], dims[2], (float)0);
	dvals[3] = GetData3(D, x1, y1, z0, dims[0], dims[1], dims[2], (float)0);
	dvals[4] = GetData3(D, x0, y0, z1, dims[0], dims[1], dims[2], (float)0);
	dvals[5] = GetData3(D, x1, y0, z1, dims[0], dims[1], dims[2], (float)0);
	dvals[6] = GetData3(D, x0, y1, z1, dims[0], dims[1], dims[2], (float)0);
	dvals[7] = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
	float dvalsX[4];
	dvalsX[0] = (1.0-r)*dvals[0] + r*dvals[1];
	dvalsX[1] = (1.0-r)*dvals[2] + r*dvals[3];
	dvalsX[2] = (1.0-r)*dvals[4] + r*dvals[5];
	dvalsX[3] = (1.0-r)*dvals[6] + r*dvals[7];
	float dvalsY[2];
	dvalsY[0] = (1.0-s)*dvalsX[0] + s*dvalsX[1];
	dvalsY[1] = (1.0-s)*dvalsX[2] + s*dvalsX[3];
	float dvalsZ;
	dvalsZ = (1.0-t)*dvalsY[0] + t*dvalsY[1];
	
	return dvalsZ;
}

float
Divergence(const vector<float>& D,
		   int x, int y, int z, 
		   const int* dims)
{
	float dval = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	float sum = 0;
	for(int j=0; j<NumNeighbors; ++j)
	{
		int x2 = x+XOffset[j];
		int y2 = y+YOffset[j];
		int z2 = z+ZOffset[j];
		float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)dval);
		sum += dval - dval2;
	}

	return sum;
}

vector<CParticle>
collectConnectedNeighborRegion(const vector<unsigned char>& L,
							   const CParticle& p,
							   const vector<float>& v,
							   float maxDistance,
							   const int* dims)
{
	vector<float> F(L.size(), (float)-1);
	vector<CParticle> vp(1, p);
	vector<CParticle> vregion(1, p);
	SetData3(F, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], (float)0);
	while(!vp.empty())
	{
		vector<CParticle> vp2;
		for(int n=0; n<vp.size(); ++n)
		{
			int x=vp[n].m_X;
			int y=vp[n].m_Y;
			int z=vp[n].m_Z;
			float fval = GetData3(F, x, y, z, dims[0], dims[1], dims[2], (float)0);
			for(int j=0; j<NumNeighbors; ++j)
			{
				int x2 = x+XOffset[j];
				int y2 = y+YOffset[j];
				int z2 = z+ZOffset[j];
				if(GetData3(L, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)0) == 0)
				{
					//outside the foreground
					continue;
				}
				float fval2 = GetData3(F, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
				if(fval2 >= 0)
				{
					//already visited.
					continue;
				}
				float dx = XOffset[j]*v[0];
				float dy = YOffset[j]*v[1];
				float dz = ZOffset[j]*v[2];
				float dd = sqrt(dx*dx + dy*dy + dz*dz);
				float fvalnew = fval + dd;
				if(fvalnew <= maxDistance)
				{
					vp2.push_back(CParticle(x2, y2, z2, fvalnew));
				}
			}
		}
		for(int i=0; i<vp2.size(); ++i)
		{
			float fval2 = GetData3(F, vp2[i].m_X, vp2[i].m_Y, vp2[i].m_Z, dims[0], dims[1], dims[2], (float)0);
			if(fval2 < 0 || fval2 > vp2[i].m_Life)
			{
				SetData3(F, vp2[i].m_X, vp2[i].m_Y, vp2[i].m_Z, dims[0], dims[1], dims[2], (float)vp2[i].m_Life);
			}
		}
		vector<CParticle> vnext;
		for(int i=0; i<vp2.size(); ++i)
		{
			float fval2 = GetData3(F, vp2[i].m_X, vp2[i].m_Y, vp2[i].m_Z, dims[0], dims[1], dims[2], (float)0);
			if(fval2 == vp2[i].m_Life)
			{
				vnext.push_back(vp2[i]);
			}
		}

		vregion.insert(vregion.end(), vnext.begin(), vnext.end());
		vp = vnext;
	}
	return vregion;
}