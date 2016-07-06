#ifndef _SCpoint_h_
#define _SCpoint_h_
//Quadratic surface patch
#ifdef MEX_DLL
#include <mex.h>
#endif

#include <iostream>
#include <fstream>
using namespace std;

#include <stdlib.h>
#include <assert.h>
#include <mytype.h>
#include <Complex.h>

/*
Define a surface patch:
f(s,t|X,Y,Z)=Z+A*s^3+B*s*t^2+C*s*t^2+D*t^3+E*s^2+F*s*t+G*t^2+H*s+P*t;
*/
class SCpoint{
public:
  SCpoint(real x=0, real y=0, real z=0, 
    real a=0, real b=0, real c=0, real d=0, real e=0,
    real f=0, real g=0, real h=0, real p=0) {
    X=x;
    Y=y;
    Z=z;
    A=a;
    B=b;
    C=c;
    D=d;
    E=e;
    F=f;
    G=g;
    H=h;
    P=p;
  }
  
  SCpoint(const SCpoint& cp) {
    X=cp.GetX();
    Y=cp.GetY();
    Z=cp.GetZ();
    A=cp.GetA();
    B=cp.GetB();
    C=cp.GetC();
    D=cp.GetD();
    E=cp.GetE();
    F=cp.GetF();
    G=cp.GetG();
    H=cp.GetH();
    P=cp.GetP();
  }
  
  const SCpoint& operator=(const SCpoint& cp) {
    X=cp.GetX();
    Y=cp.GetY();
    Z=cp.GetZ();
    A=cp.GetA();
    B=cp.GetB();
    C=cp.GetC();
    D=cp.GetD();
    E=cp.GetE();
    F=cp.GetF();
    G=cp.GetG();
    H=cp.GetH();
    P=cp.GetP();
    
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
  inline real GetF() const {return F;}
  inline real GetG() const {return G;}
  inline real GetH() const {return H;}
  inline real GetP() const {return P;}
  inline void SetX(real x) {X=x;}
  inline void SetY(real y) {Y=y;}
  inline void SetZ(real z) {Z=z;}
  inline void SetA(real a) {A=a;}
  inline void SetB(real b) {B=b;}
  inline void SetC(real c) {C=c;}
  inline void SetD(real d) {D=d;}
  inline void SetE(real e) {E=e;}
  inline void SetF(real f) {F=f;}
  inline void SetG(real g) {G=g;}
  inline void SetH(real h) {H=h;}
  inline void SetP(real p) {P=p;}
  
  real Evaluate(real s, real t) const {
    return Z+A*s*s*s+B*s*s*t+C*s*t*t+D*t*t*t+E*s*s+F*s*t+G*t*t+H*s+P*t;
  }
  
  Complex Tangent(real s, real t) const {
    return Complex(3.*A*s*s+2*B*s*t+C*t*t+2*E*s+F*t+H,\
       B*s*s+2*C*s*t+3*D*t*t+F*s+2*G*t+P);
  }

  CCurvature Curvature(real s, real t);
  
  // Artithmetic overloaded operators
  /*SCpoint operator+(const SCpoint& s) const {
    SCpoint p(X+s.GetX(),Y+s.GetY(),Z+s.GetZ(),\
      A+s.GetA(),B+s.GetB(),C+s.GetC(),D+s.GetD(),E+s.GetE(),\
      F+s.GetF(),G+s.GetG(),H+s.GetH(),P+s.GetP());
    return p;
  }
  SCpoint operator-(const SCpoint& s) const {
    SCpoint p(X-s.GetX(),Y-s.GetY(),Z-s.GetZ(),\
      A-s.GetA(),B-s.GetB(),C-s.GetC(),D-s.GetD(),E-s.GetE(),\
      F-s.GetF(),G-s.GetG(),H-s.GetH(),P-s.GetP());
    return p;
  }
  
  void operator+=(const SCpoint& s) {
    X+=s.GetX(); Y+=s.GetY(); Z+=s.GetZ(); 
    A+=s.GetA(); B+=s.GetB(); C+=s.GetC(); D+=s.GetD(); E+=s.GetE();
    F+=s.GetF(); G+=s.GetG(); H+=s.GetH(); P+=s.GetP();
  }
  void operator-=(const SCpoint& s) {
    X-=s.GetX(); Y-=s.GetY(); Z-=s.GetZ(); 
    A-=s.GetA(); B-=s.GetB(); C-=s.GetC(); D-=s.GetD(); E-=s.GetE();
    F-=s.GetF(); G-=s.GetG(); H-=s.GetH(); P-=s.GetP();
  }
  
  //real-complex operators
  SCpoint operator*(real v) const {
    return SCpoint(v*X,v*Y,v*Z,v*A,v*B,v*C,v*D,v*E,v*F,v*G,v*H,v*P);
  }
  void operator*=(real v) {
    X*=v; Y*=v; Z*=v; A*=v; B*=v; C*=v; D*=v; E*=v; F*=v; G*=v; H*=v; P*=v;
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
   real F;
   real G;
   real H;
   real P;
};

// Standard Output
ostream& operator<<(ostream& s, const SCpoint& v);

istream& operator>>(istream& s, SCpoint& v);

typedef myImage<SCpoint> SCpointImage;
const int MatDimC=10;

SCpointImage
InitializeSurfaceC(const RealImage& image);

int
ReadGammaMatricesC(real* g1,real* g2,real* g3,real* g4, 
                  real& tau, real& lambda,
                  char* filename);

void
UpdateCubicSurface(SCpointImage& surf, 
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda);

void
UpdateCubicSurfaceBreak(SCpointImage& surf,
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda, real thres);

void
UpdateCubicSurfaceMask(SCpointImage& surf,
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda, 
						const ByteImage& mask);

#endif /* SCpoint_h */
