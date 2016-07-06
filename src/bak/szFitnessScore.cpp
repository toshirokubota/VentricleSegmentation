
//#include <mex.h>
#include <string.h>
#include <math.h>
#include <algorithm>
using namespace std;
#include <stdio.h>

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
//#include <myNeighborOp.h>
//#include <DistanceTransform.h>
//#include <LocalExtrema.h>
//#include <ConnectedComponent.h>
#include <szMiscOperations.h>

const int NumFeatures = 13;
const int NumFeaturesNoConvexHull = 12;

int
NumberOfFeatures()
{
	return NumFeatures;
}

/*
Compute the diameter of the segmentation in the current axial slice.
*/
int
axialDiameterEllipse(double& major,            //OUTPUT - major diameter
                     double& minor,            //OUTPUT - minor diameter
                     vector<unsigned char>& C, //segmented volume
                     int slice,                //slice number
                     vector<float>& v,        //voxel size
                     int ndim,
                     const int* dims) 
{

  int numpixels = dims[0]*dims[1];
  int offset = slice*numpixels;

  vector<int> vx;
  vector<int> vy;
  //collect non-zero voxels that are on the boundary of the 2D component
  int i,j,k;
  for(i=0, k=offset; i<dims[1]; ++i) 
  {
    for(j=0; j<dims[0]; ++j, k++) 
	{
      if(C[k]) 
	  {
        vx.push_back(j);
        vy.push_back(i);
      }
    }
  }
  int ret = ComputeAxialDiameterEllipse(major,minor,vx,vy,v);
    
  return ret;
}

int
axialDiameterDirect(double& major,            //OUTPUT - major diameter
                    double& minor,            //OUTPUT - minor diameter
                    vector<unsigned char>& C, //segmented volume
                    int slice,                //slice number
                    vector<float>& v,        //voxel size
                    int ndim,
                    const int* dims) 
{

  int numpixels = dims[0]*dims[1];
  int offset = slice*numpixels;

  vector<int> vx;
  vector<int> vy;
  //collect non-zero voxels that are on the boundary of the 2D component
  int i,j,k;
  for(i=0, k=offset; i<dims[1]; ++i) 
  {
    for(j=0; j<dims[0]; ++j, k++) 
	{
      if(C[k]) 
	  {
        if(onSurfaceOnSlice(C,j,i,slice,dims)) 
		{
          vx.push_back(j);
          vy.push_back(i);
        }
      }
    }
  }
  int ret = ComputeAxialDiameterDirect(major,minor,vx,vy,v);
    
  return ret;
}

/*
vector<double>
fitnessScoreNonSolitary(vector<unsigned char>& nC, //non-convex region - before convex Hull
                        vector<unsigned char>& C,  //convex region - after convex Hull
                        vector<unsigned char>& L,  //original foreground segmentation
                        vector<float>& D,         //distance map
                        vector<float>& Sp,        //special sphericity measure
                        vector<float>& v,         //voxel size
                        int ndim,
                        const int* dims) 
{

  //printf("Enter fitnessScoreNonSolitary()\n");
  if(v.size()<3) 
  {
    printf("Error: There are not enough voxel size elements...");
    return vector<double>(NumFeatures,0);
  }
  int i;
  int nvoxels=numberOfElements(ndim,dims);

  double maxSp=0;
  double maxD=0;
  int vols=0; //# voxels in non-convex (nC) volume
  int difvols=0; //# voxels in convex (C) volume but not in non-convex
  int fitCount=0; //# voxels that are on surface of nC and on surface of L
  int unfitCount=0; //# voxels that are on surface of nC but not on surface of L
  int maxIdx=0;
  for(i=0; i<nvoxels; ++i) 
  {
    if(nC[i]) {
      vols++;
      if(maxSp<Sp[i]) 
	  {
        maxSp=Sp[i];
        maxIdx=i;
      }
      if(maxD<D[i])
        maxD=D[i];
      if(onSurface3(nC,i,dims)) 
	  {
        if(onSurface3(L,i,dims))
          fitCount++;
        else
          unfitCount++;
      }
    }
    else if(C[i]) 
	{
      difvols++;
    }
  }
  //compute sphericity with a larger neighborhood
  int ix,iy,iz;
  Ind2Sub(ix,iy,iz,maxIdx,dims);
  float spval=EvaluateSphericityValue(D,ix,iy,iz,2,2,2,dims);


  //printf("Computing axial diameter.\n");
  double major1=0, minor1=0, major2=0, minor2=0;
  for(i=0; i<dims[2]; ++i) 
  {
    double d1=0, d2=0, d3=0, d4=0;
    int bRet = axialDiameterDirect(d1,d2,nC,i,v,ndim,dims);
    if(bRet>0 && d1*d2>major1*minor1) 
	{
      major1 = d1;
      minor1 = d2;
    }
    int bRet2 = axialDiameterEllipse(d3,d4,nC,i,v,ndim,dims);
    if(bRet2>0 && d3*d4>major2*minor2) 
	{
      major2 = d3;
      minor2 = d4;
    }
  }
  //mexPrintf("Done computing axial diameter. %d %f %f\n", vols, major, minor);
  if(vols==0)
    return vector<double>(NumFeatures,0);

  vector<double> score;
  score.push_back(vols);
  score.push_back(vols*v[0]*v[1]*v[2]);
  score.push_back(maxSp);
  score.push_back(spval);
  score.push_back(maxD);
  score.push_back((double)difvols/vols);
  if(fitCount || unfitCount) 
  {
    score.push_back((double)fitCount/(fitCount+unfitCount));
    //score.push_back((double)unfitCount/(fitCount+unfitCount));
  }
  else 
  { //this should not happen when vols>0
    score.push_back(0);
    //score.push_back(0);
  }
  score.push_back((major1+minor1)/2);
  score.push_back(sqrt(major1*minor1));
  if(major1>0)
    score.push_back(minor1/major1);
  else
    score.push_back(0);
  score.push_back((major2+minor2)/2);
  score.push_back(sqrt(major2*minor2));
  if(major2>0)
    score.push_back(minor2/major2);
  else
    score.push_back(0);

  //printf("Exit fitnessScoreNonSolitary()\n");
  return score;
}
*/

/*
vector<double>
fitnessScoreNonSolitaryNoConvexHull(vector<unsigned char>& nC, //non-convex region - before convex Hull
                                    vector<unsigned char>& L,  //original foreground segmentation
                                    vector<float>& D,         //distance map
                                    vector<float>& Sp,        //special sphericity measure
                                    vector<float>& v,         //voxel size
                                    int ndim,
                                    const int* dims) 
{

  //printf("Enter fitnessScoreNonSolitary()\n");
  if(v.size()<3) 
  {
    printf("Error: There are not enough voxel size elements...");
    return vector<double>(NumFeaturesNoConvexHull,0);
  }
  int i;
  int nvoxels=numberOfElements(ndim,dims);

  double maxSp=0;
  double maxD=0;
  int vols=0; //# voxels in non-convex (nC) volume
  int difvols=0; //# voxels in convex (C) volume but not in non-convex
  int fitCount=0; //# voxels that are on surface of nC and on surface of L
  int unfitCount=0; //# voxels that are on surface of nC but not on surface of L
  int maxIdx=0;
  for(i=0; i<nvoxels; ++i) 
  {
    if(nC[i]) {
      vols++;
      if(maxSp<Sp[i]) 
	  {
        maxSp=Sp[i];
        maxIdx=i;
      }
      if(maxD<D[i])
        maxD=D[i];
      if(onSurface3(nC,i,dims)) 
	  {
        if(onSurface3(L,i,dims))
          fitCount++;
        else
          unfitCount++;
      }
    }
  }
  //compute sphericity with a larger neighborhood
  int ix,iy,iz;
  Ind2Sub(ix,iy,iz,maxIdx,dims);
  float spval=EvaluateSphericityValue(D,ix,iy,iz,2,2,2,dims);

  //printf("Computing axial diameter.\n");
  double major1=0, minor1=0, major2=0, minor2=0;
  for(i=0; i<dims[2]; ++i) 
  {
    double d1=0, d2=0, d3=0, d4=0;
    int bRet = axialDiameterDirect(d1,d2,nC,i,v,ndim,dims);
    if(bRet>0 && d1*d2>major1*minor1) 
	{
      major1 = d1;
      minor1 = d2;
    }
    int bRet2 = axialDiameterEllipse(d3,d4,nC,i,v,ndim,dims);
    if(bRet2>0 && d3*d4>major2*minor2) 
	{
      major2 = d3;
      minor2 = d4;
    }
  }
  //mexPrintf("Done computing axial diameter. %d %f %f\n", vols, major, minor);
  if(vols==0)
    return vector<double>(NumFeatures,0);

  vector<double> score;
  score.push_back(vols);
  score.push_back(vols*v[0]*v[1]*v[2]);
  score.push_back(maxSp);
  score.push_back(spval);
  score.push_back(maxD);
  if(fitCount || unfitCount) 
  {
    score.push_back((double)fitCount/(fitCount+unfitCount));
    //score.push_back((double)unfitCount/(fitCount+unfitCount));
  }
  else 
  { //this should not happen when vols>0
    score.push_back(0);
    //score.push_back(0);
  }
  score.push_back((major1+minor1)/2);
  score.push_back(sqrt(major1*minor1));
  if(major1>0)
    score.push_back(minor1/major1);
  else
    score.push_back(0);
  score.push_back((major2+minor2)/2);
  score.push_back(sqrt(major2*minor2));
  if(major2>0)
    score.push_back(minor2/major2);
  else
    score.push_back(0);

  //printf("Exit fitnessScoreNonSolitary()\n");
  return score;
}
*/

/*
vector<double>
fitnessScoreDegenerate() 
{

  vector<double> score(NumFeatures,0);

  return score;
}
*/

/*
vector<double>
fitnessScoreDegenerateNoConvexHull() 
{

  vector<double> score(NumFeaturesNoConvexHull,0);

  return score;
}
*/

vector<double>
fitnessScore(vector<unsigned char>& B, //non-convex region - before convex Hull
			 vector<float>& D,         //distance map
			 vector<float>& Sp,        //special sphericity measure
			 vector<float>& v,         //voxel size
			 int ndim,
			 const int* dims) 
{

  //printf("Enter fitnessScoreNonSolitary()\n");
  if(v.size()<3) 
  {
    printf("Error: There are not enough voxel size elements...");
    return vector<double>(NumFeaturesNoConvexHull,0);
  }
  int i;
  int nvoxels=numberOfElements(ndim,dims);

  double maxSp=0;
  double maxD=0;
  int vols=0; //# voxels in non-convex (nC) volume
  int maxDIdx = 0, maxSpIdx = 0;
  for(i=0; i<nvoxels; ++i) 
  {
    if(B[i]) {
      vols++;
      if(maxSp<Sp[i]) 
	  {
        maxSp=Sp[i];
        maxSpIdx=i;
      }
      if(maxD<D[i])
	  {
        maxD=D[i];
		maxDIdx = i;
	  }
    }
  }
  if(vols==0)
    return vector<double>(NumFeatures,0);

  //compute sphericity with a larger neighborhood
  float maxSphericity = MaximumSphericityValue(2, 2, 2);
  int ix,iy,iz;
  Ind2Sub(ix,iy,iz,maxSpIdx,dims);
  float spval1=EvaluateSphericityValue(D,ix,iy,iz,2,2,2,dims) / maxSphericity;
  Ind2Sub(ix,iy,iz,maxDIdx,dims);
  float spval2=EvaluateSphericityValue(D,ix,iy,iz,2,2,2,dims) / maxSphericity;


  //go through slice by slice and store the maximum diameter
  double major1=0, minor1=0, major2=0, minor2=0;
  for(i=0; i<dims[2]; ++i) {
    double d1=0, d2=0, d3=0, d4=0;
    int bRet = axialDiameterDirect(d1,d2,B,i,v,ndim,dims);
    if(bRet>0 && d1*d2>major1*minor1) {
      major1 = d1;
      minor1 = d2;
    }
    int bRet2 = axialDiameterEllipse(d3,d4,B,i,v,ndim,dims);
    if(bRet2>0 && d3*d4>major2*minor2) {
      major2 = d3;
      minor2 = d4;
    }
  }
  //compute the dimater at the maximum distant point
  double major3=0, minor3=0, major4=0, minor4=0;
  Ind2Sub(ix,iy,iz,maxDIdx,dims);
  {
	  double d1=0, d2=0, d3=0, d4=0;
	  int bRet = axialDiameterDirect(d1,d2,B,iz,v,ndim,dims);
	  if(bRet>0) {
		  major3 = d1;
		  minor3 = d2;
	  }
	  int bRet2 = axialDiameterEllipse(d3,d4,B,iz,v,ndim,dims);
	  if(bRet2>0) {
		  major4 = d3;
		  minor4 = d4;
	  }
  }

  vector<double> score;
  score.push_back(vols);
  score.push_back(vols*v[0]*v[1]*v[2]);
  score.push_back(spval1);
  score.push_back(spval2);
  score.push_back(maxD);
  score.push_back((major1+minor1)/2);
  score.push_back(major1);
  score.push_back((major2+minor2)/2);
  score.push_back(major2);
  score.push_back((major3+minor3)/2);
  score.push_back(major3);
  score.push_back((major4+minor4)/2);
  score.push_back(major4);

  //printf("Exit fitnessScoreNonSolitary()\n");
  return score;
}

