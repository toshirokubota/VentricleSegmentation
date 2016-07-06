#ifndef _DQpoint_h_
#define _DQpoint_h_
//Quadratic density patch

#include <iostream>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <assert.h>
#include "mytype.h"
#include "Complex.h"

/*
Define a surface patch:
f(s,t|X,Y,Z)=Z+A*s^2+B*s+C*t^2+D*t+E*s*t;
*/
class DQpoint{
public:
  DQpoint(real x=0, real y=0, real z=0, 
    real a=0, real b=0, real c=0, real d=0, real e=0, real f=0) {
    X=x;
	Y=y;
	Z=z;
    A=a;
    B=b;
    C=c;
    D=d;
    E=e;
	F=f;
  }
  
  DQpoint(const DQpoint& cp) {
    X=cp.GetX();
    Y=cp.GetY();
    Z=cp.GetZ();
    A=cp.GetA();
    B=cp.GetB();
    C=cp.GetC();
    D=cp.GetD();
    E=cp.GetE();
    F=cp.GetF();
  }
  
  const DQpoint& operator=(const DQpoint& cp) {
    X=cp.GetX();
    Y=cp.GetY();
    Z=cp.GetZ();
    A=cp.GetA();
    B=cp.GetB();
    C=cp.GetC();
    D=cp.GetD();
    E=cp.GetE();
    F=cp.GetF();
    
    return *this;
  }
  
  inline real GetX() const {return X;}
  inline real GetY() const {return Y;}
  inline real GetZ() const {return Z;}
  inline real GetA() const {return A;}
  inline real GetB() const {return B;}
  inline real GetC() const {return C;}
  inline real GetD() const {return D;}
  inline real GetE() const {return E;}
  inline real GetF() const {return E;}
  inline void SetX(real x) {
    X=x;
  }
  inline void SetY(real y) {
    Y=y;
  }
  inline void SetZ(real z) {
    Z=z;
  }
  inline void SetA(real a) {
    A=a;
  }
  inline void SetB(real b) {
    B=b;
  }
  inline void SetC(real c) {
    C=c;
  }
  inline void SetD(real d) {
    D=d;
  }
  inline void SetE(real e) {
    E=e;
  }
  inline void SetF(real f) {
    F=f;
  }
  
 protected:
   real X;
   real Y;
   real Z;
   real A;
   real B;
   real C;
   real D;
   real E;
   real F;
};

ostream& operator<<(ostream& s, const DQpoint& v);

istream& operator>>(istream& s, DQpoint& v);

typedef myImage<DQpoint> DQpointImage;
const int MatDimDQ=7;

void
InitializeDensityQ(DQpointImage& res, const RealImage& image);

void
InitializeDensityQ(DQpointImage& res, const vector<real>& image, const int* dims);

void
UpdateQuadraticDensity(DQpointImage& density, 
					   real r, real t);

void
UpdateQuadraticDensityMask(DQpointImage& density, 
						   real r, real t,
						   const vector<unsigned char>& M);

#endif /* DQpoint_h */
