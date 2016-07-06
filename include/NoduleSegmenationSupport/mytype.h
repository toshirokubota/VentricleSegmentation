#ifndef ___MYTYPE_H___
#define ___MYTYPE_H___
#include<mex.h>

#include <vector>
using namespace std;
#include <myImage.h>

typedef float real;
typedef unsigned char uchar;

typedef vector<real> vReal;
typedef vector<int> vInt;

typedef myImage<unsigned char> ByteImage;
typedef myImage<real> RealImage;
typedef myImage<int> IntImage;

class CGradient
{
public:
	real x;
	real y;
};

class CCurvature
{
public:
	real H;
	real K;
	real k1;
	real k2;
};

#endif /* ___MYTYPE_H___ */
