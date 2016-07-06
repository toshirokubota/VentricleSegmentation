#ifndef ___MY_EIGEN_H___
#define ___MY_EIGEN_H___

#include <vector>
using namespace std;

bool
Eigen(const vector<double>& C, 
	  vector<double>& eigenVal, 
	  vector<double>& eigenVec, 
	  int ndim);

vector<double>
computeEigenValues(const vector<double>& M, int ndim);

#endif /* ___MY_EIGEN_H___ */