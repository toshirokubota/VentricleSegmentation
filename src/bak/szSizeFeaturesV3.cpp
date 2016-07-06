
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
#include <szEigen.h>
#include <szFeatureClass.h>

inline int
NumSizeFeatures()
{
	const int NumFeatures = 7;
	return NumFeatures;
}

/*
Compute the diameter of the segmentation in the current axial slice.
*/
int
axialDiameterEllipse3D(double& diameter1,            //OUTPUT - major diameter
					   double& diameter2,            //OUTPUT - minor diameter
					   double& diameter3,            //OUTPUT - minor diameter
					   const vector<unsigned char>& C, //segmented volume
					   const vector<float>& v,        //voxel size
					   int ndim,
					   const int* dims) 
{
	vector<double> vx;
	vector<double> vy;
	vector<double> vz;
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(C, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					vx.push_back(v[0]*k);
					vy.push_back(v[1]*j);
					vz.push_back(v[2]*i);
				}
			}
		}
	}
	vector<double> vcov;
	diameter1 = diameter2 = diameter3 = 0;
	if(vx.size()>0)
	{
		vcov = computeCovariance(vx, vy, vz);
		vector<double> veigen = computeEigenValues(vcov, ndim);
		if(veigen.size() == ndim)
		{
			diameter1 = veigen[2];
			diameter2 = veigen[1];
			diameter3 = veigen[0];
		}
		printf("axialDiameterEllipse3D: diameters = %f, %f, %f\n", 
			diameter1, diameter2, diameter3);
		return 1;
	}
	else
	{
		printf("axialDiameterEllipse3D: Failed.\n", 
			diameter1, diameter2, diameter3);
		return 0;
	}
}

/*
Compute the diameter of the segmentation in the current axial slice.
*/
int
axialDiameterEllipse(double& major,            //OUTPUT - major diameter
                     double& minor,            //OUTPUT - minor diameter
                     const vector<unsigned char>& C, //segmented volume
                     int slice,                //slice number
                     const vector<float>& v,        //voxel size
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
                    const vector<unsigned char>& C, //segmented volume
                    int slice,                //slice number
                    const vector<float>& v,        //voxel size
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

vector<CFeature>
computeSizeFeaturesV3(const vector<unsigned char>& B, //nodule segmentation
					  const vector<float>& v,         //voxel size
					  int ndim,
					  const int* dims) 
{
  int i;
  int nvoxels=numberOfElements(ndim,dims);

  int vols=0; //# voxels in non-convex (nC) volume
  for(i=0; i<nvoxels; ++i) 
  {
    if(B[i]) {
      vols++;
    }
  }

  vector<CFeature> features(NumSizeFeatures());
  if(vols==0)
    return features;

  //Diameter features
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

  features[0].Value = vols * v[0] * v[1] * v[2];
  features[0].Name = "SegmentationVolume";
  features[1].Value = (major1+minor1)/2;
  features[1].Name = "ELCAPDiameterMean";
  if(major1>0)
	features[2].Value = minor1/major1;
  features[2].Name = "ELCAPDiameterRatio";
  features[3].Value = (major2+minor2)/2;
  features[3].Name = "EllipseMeanRadius";
  if(major2>0)
	features[4].Value = minor2/major2;
  features[4].Name = "EllipseRadiiRatio";

  double a1, a2, a3;
  axialDiameterEllipse3D(a1, a2, a3, B, v, ndim, dims);
  if(a1)
  {
	  features[5].Value = a3/a1;
  }
  features[5].Name = "EllipsoidRadiiRatio1";
  if(a2)
  {
	  features[6].Value = a3/a2;
  }
  features[6].Name = "EllipsoidRadiiRatio2";

  return features;
}
