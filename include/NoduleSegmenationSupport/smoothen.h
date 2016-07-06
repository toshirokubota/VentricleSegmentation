#ifndef ___SMOOTHEN_H___
#define ___SMOOTHEN_H___
#ifdef MEX_DLL
#include <mex.h>
#endif

#include <vector>
using namespace std;

#include <mytype.h>
#include <SCpoint.h>
//typedef float real;

void
smoothen3(vector<float>& A, 
		  const vector<unsigned char>& M,
		 int niter,
		 real* gn, real* gs, real* gw, real* ge,
		 real tau, real lambda,
		 const int* dims);

void
smoothen3D(vector<float>& A, 
		  const vector<unsigned char>& M,
		 int niter,
		 real tau, real lambda,
		 const int* dims);

void
smoothen3L(vector<float>& A, 
		   const vector<unsigned char>& M,
		   int niter,
		   float alpha,
		   const int* dims);

#endif /* ___SMOOTHEN_H___ */