#ifndef _Complex_h_
#define _Complex_h_

//#include <iostream>
//using namespace std;

#include <math.h>
#include "mytype.h"

class Complex {
 public:
  Complex(real r=0, real img=0) {Real=r; Imag=img;}
  Complex(const Complex& comp_num);

  const Complex& operator=(const Complex& comp_num);

  real GetReal() const {return Real;}
  real GetImag() const {return Imag;}
  void SetReal(const real r) {Real = r;}
  void SetImag(const real im) {Imag = im;}

  Complex operator+(const Complex& comp_num) const;
  Complex operator-(const Complex& comp_num) const;
  Complex operator*(const Complex& comp_num) const;

  void operator+=(const Complex& comp_num);
  void operator-=(const Complex& comp_num);
  void operator*=(const Complex& comp_num);

  //real-complex operators
  Complex operator*(real real_num) const;
  void operator*=(real real_num);

  real Magnitude() const;
  real Power() const;
  real Phase() const;
  real PhaseB() const;

 private:
  real Real;
  real Imag;
};

enum Complex2Real_Op {Complex2Real_Real, Complex2Real_Imag,
                      Complex2Real_Mag, Complex2Real_Phase};

istream& operator>>(istream& s, Complex& comp_num);
ostream& operator<<(ostream& s, const Complex& comp_num);

Complex operator* (real real_num, const Complex& comp_num);

typedef myImage<Complex> ComplexImage;

#endif /* Complex_h */
