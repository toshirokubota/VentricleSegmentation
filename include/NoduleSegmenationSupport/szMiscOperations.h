#ifndef _SZ_MISC_OPERATIONS_H_
#define _SZ_MISC_OPERATIONS_H_
#include <vector>
using namespace std;
#include <szParticle.h>

int GetOffsetX6(int n);

int GetOffsetY6(int n);

int GetOffsetZ6(int n);

int GetOffsetX10(int n);

int GetOffsetY10(int n);

int GetOffsetZ10(int n);

int GetOffsetX26(int n);

int GetOffsetY26(int n);

int GetOffsetZ26(int n);

#include <szIndexedData.h>

vector<int>
indexedSort(vector<float>& vvalues);

bool
onSurface3(const vector<unsigned char>& A, 
           int x, int y, int z, 
           const int* dims);
bool
onSurface3(const vector<int>& A, 
           int x, int y, int z, 
           const int* dims);
bool
onBoundary3(const vector<unsigned char>& A, 
			int x, int y, int z, 
			const int* dims);

bool
onSurface3(const vector<unsigned char>& A, int id, const int* dims);

bool
onAdjacentSurface3(const vector<unsigned char>& A, 
                   int x, int y, int z, 
                   const int* dims);

bool
onSurfaceOnSlice(const vector<unsigned char>& A, 
                 int x, int y, int z, 
                 const int* dims);

int
ComputeAxialDiameterDirect(double& major,       //OUTPUT - major diameter
                           double& minor,       //OUTPUT - minor diameter
                           const vector<int>& vx,     //x coordinates
                           const vector<int>& vy,     //y cordinates
                           const vector<float>& v);   //voxel size

int
ComputeAxialDiameterEllipse(double& major,       //OUTPUT - major diameter
                            double& minor,       //OUTPUT - minor diameter
                            const vector<int>& vx,     //x coordinates
                            const vector<int>& vy,     //y cordinates
                            const vector<float>& v);   //voxel size

float
EvaluateSphericityValue(const vector<float>& D,
                        int x,
                        int y, 
                        int z,
                        int widthX,
                        int widthY,
                        int widthZ,
                        const int* dims);

float
MaximumSphericityValue(int widthX,
                       int widthY,
                       int widthZ);

void
selectCluster(vector<unsigned char>& L, 
              int x, int y, int z, 
              const int* dims);


vector<int>
representativePosition(const vector<int>& pos);

vector<unsigned char>
refineSegmentation(vector<unsigned char>& S, //segmentation
				   const vector<unsigned char>& L, //foreground
				   const vector<int>& cn, //a causal neighborhood as a structural element
				   int ndim, 
				   const int* dims //dimension of the volume
				   );
/*
compare two foreground segmentations to pick the better one.
*/
int 
compareForegroundFitness(const vector<float>& vSp1,
						 const vector<int>& vseeds1,
						 const vector<float>& vSp2,
						 const vector<int>& vseeds2,
						 const int* dims);

int 
compareForegroundFitness2(const vector<float>& vSp1,
						  int noduleType1,
						  const vector<float>& vSp2,
						  int noduleType2);

vector<double>
computeCovariance(const vector<double>& vx,
				  const vector<double>& vy,
				  const vector<double>& vz);

vector<double>
computeCentroid(const vector<unsigned char>& S, //segmentation
				int ndim, //should be 3
				const int* dims);

void
FillBackgroundHoles(vector<unsigned char>& B,
					int minSize,
					int label,
					const int* dims);

float
trilinearInterpolation(const vector<float>& D,
					   float x, 
					   float y, 
					   float z,
					   const int* dims);

float
Divergence(const vector<float>& D,
		   int x, int y, int z, 
		   const int* dims);

vector<CParticle>
collectConnectedNeighborRegion(const vector<unsigned char>& L,
							   const CParticle& p,
							   const vector<float>& v,
							   float maxDistance,
							   const int* dims);

#endif /* _SZ_MISC_OPERATIONS_H_ */