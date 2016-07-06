#ifndef _SQpoint_h_
#define _SQpoint_h_
//Quadratic surface patch

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
class SQpoint
{
public:
  SQpoint(real x=0, real y=0, real z=0, 
    real a=0, real b=0, real c=0, real d=0, real e=0) 
  {
    X=x;
    Y=y;
    Z=z;
    A=a;
    B=b;
    C=c;
    D=d;
    E=e;
  }
  
  SQpoint(const SQpoint& cp) 
  {
    X=cp.GetX();
    Y=cp.GetY();
    Z=cp.GetZ();
    A=cp.GetA();
    B=cp.GetB();
    C=cp.GetC();
    D=cp.GetD();
    E=cp.GetE();
  }
  
  const SQpoint& operator=(const SQpoint& cp) 
  {
    X=cp.GetX();
    Y=cp.GetY();
    Z=cp.GetZ();
    A=cp.GetA();
    B=cp.GetB();
    C=cp.GetC();
    D=cp.GetD();
    E=cp.GetE();
    
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
  inline void SetX(real x) 
  {
    X=x;
  }
  inline void SetY(real y) 
  {
    Y=y;
  }
  inline void SetZ(real z) 
  {
    Z=z;
  }
  inline void SetA(real a) 
  {
    A=a;
  }
  inline void SetB(real b) 
  {
    B=b;
  }
  inline void SetC(real c) 
  {
    C=c;
  }
  inline void SetD(real d) 
  {
    D=d;
  }
  inline void SetE(real e) 
  {
    E=e;
  }
  
  /*real Evaluate(real s, real t) const 
  {
    return Z+A*s*s+B*s+C*t*t+D*t+E*s*t;
  }*/
  
  Complex Tangent(real s, real t) const 
  {
    return Complex(2.*A*s+B+E*t,2.*C*t+D+E*s);
  }
  
  CCurvature Curvature(real s, real t) const;

  // Artithmetic overloaded operators
  /*SQpoint operator+(const SQpoint& s) const 
  {
    SQpoint p(X+s.GetX(),Y+s.GetY(),Z+s.GetZ(),\
      A+s.GetA(),B+s.GetB(),C+s.GetC(),D+s.GetD(),E+s.GetE());
    return p;
  }
  SQpoint operator-(const SQpoint& s) const 
  {
    SQpoint p(X-s.GetX(),Y-s.GetY(),Z-s.GetZ(),\
      A-s.GetA(),B-s.GetB(),C-s.GetC(),D-s.GetD(),E-s.GetE());
    return p;
  }
  void operator+=(const SQpoint& s) 
  {
    X+=s.GetX(); Y+=s.GetY(); Z+=s.GetZ(); 
    A+=s.GetA(); B+=s.GetB(); C+=s.GetC(); D+=s.GetD(); E+=s.GetE();
  }
  void operator-=(const SQpoint& s) 
  {
    X-=s.GetX(); Y-=s.GetY(); Z-=s.GetZ(); 
    A-=s.GetA(); B-=s.GetB(); C-=s.GetC(); D-=s.GetD(); E-=s.GetE();
  }
  
  //real-complex operators
  SQpoint operator*(real v) const 
  {
    return SQpoint(v*X,v*Y,v*Z,v*A,v*B,v*C,v*D,v*E);
  }
  void operator*=(real v) 
  {
    X*=v; Y*=v; Z*=v; A*=v; B*=v; C*=v; D*=v; E*=v;
  }*/
  
 protected:
   real X;
   real Y;
   real Z;
   real A;
   real B;
   real C;
   real D;
   real E;
};

ostream& operator<<(ostream& s, const SQpoint& v);

istream& operator>>(istream& s, SQpoint& v);

typedef myImage<SQpoint> SQpointImage;
const int MatDimQ=6;

SQpointImage
InitializeSurfaceQ(const RealImage& image);

int
ReadGammaMatricesQ(real* g1,real* g2,real* g3,real* g4, 
                  real& tau, real& lambda,
                  char* filename);

void
UpdateQuadraticSurface(SQpointImage& surf, 
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda);

void
UpdateQuadraticSurface(SQpointImage& surf, 
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda,
						const ByteImage& mask);

void
UpdateQuadraticSurfaceBreak(SQpointImage& surf,
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda, real thres);

#endif /* SQpoint_h */
