#ifndef ___DISTANCE_TRANSFORM_H___
#define ___DISTANCE_TRANSFORM_H___

#include <string.h>
#include <math.h>
#include <vector>
using namespace std;

#include "szMexUtility.h"
#include "szMexUtilityTemplate.h"

const int DistanceTransformMode_Euclid = 1;
const int DistanceTransformMode_QuasiEuclid = 2;
const int DistanceTransformMode_CityBlockEuclid = 3;
const int DistanceTransformMode_CheckerBoardEuclid = 4;

void
DistanceTransform(vector<double>& D, 
                  const vector<unsigned char>& L, 
                  int mode,
                  int ndim, 
                  const int* dims);

void
DistanceTransformF(vector<float>& Df, 
                  const vector<unsigned char>& L, 
                  int mode,
                  int ndim, 
                  const int* dims);

/*
Invert the foreground and the background before computing the distance map
*/
void
DistanceTransformInv(vector<double>& D, 
					 const vector<unsigned char>& L, 
					 int mode,
					 int ndim, 
					 const int* dims);

/*
Invert the foreground and the background before computing the distance map
*/
void
DistanceTransformFInv(vector<float>& Df, 
					  const vector<unsigned char>& L, 
					  int mode,
					  int ndim, 
					  const int* dims);

#endif /* ___DISTANCE_TRANSFORM_H___ */