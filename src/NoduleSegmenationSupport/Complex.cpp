//#include "stdafx.h"
#include <Complex.h>

//
// Copy constructor
//
Complex::Complex(const Complex& comp_num) {

  Real = comp_num.GetReal();
  Imag = comp_num.GetImag();
}

//
// Assignment operator
//
const Complex& Complex::operator=(const Complex& comp_num) {

  if(this != &comp_num) {
    Real = comp_num.GetReal();
    Imag = comp_num.GetImag();
  }
  return *this;
}

Complex Complex::operator-(const Complex& comp_num) const {
  Complex res;

  res.SetReal(Real - comp_num.GetReal());
  res.SetImag(Imag - comp_num.GetImag());

  return res;
}

Complex Complex:: operator+(const Complex& comp_num) const {
  Complex res;

  res.SetReal(Real + comp_num.GetReal());
  res.SetImag(Imag + comp_num.GetImag());

  return res;
}

Complex Complex:: operator *(const Complex& comp_num) const {
  Complex res;

  res.SetReal(Real*comp_num.GetReal() - Imag*comp_num.GetImag());
  res.SetImag(Real*comp_num.GetImag() + Imag*comp_num.GetReal());

  return res;
}

void Complex:: operator+=(const Complex& cmpl) {

  Real += cmpl.GetReal();
  Imag += cmpl.GetImag();
}

void Complex:: operator-=(const Complex& cmpl) {

  Real -= cmpl.GetReal();
  Imag -= cmpl.GetImag();
}

void Complex:: operator*=(const Complex& cmpl) {
  Complex res;
  real tmp_real;

  tmp_real = Real*cmpl.GetReal() - Imag*cmpl.GetImag();
  Imag = Real*cmpl.GetImag() + Imag*cmpl.GetReal();
  Real = tmp_real;
}

//
// real-complex operators.
// currently, only the multiplication operator (*) is supported.
// do others make sense to have?
//
Complex Complex:: operator *(real real_num) const {
  Complex res;

  res.SetReal(Real*real_num);
  res.SetImag(Imag*real_num);

  return res;
}

void Complex:: operator*=(real real_num) {
  Real *= real_num;
  Imag *= real_num;
}

Complex operator* (real real_num, const Complex& comp_num) {
  Complex result;
  result.SetReal(real_num*comp_num.GetReal());
  result.SetImag(real_num*comp_num.GetImag());

  return result;
}

//
// Complex to real operations
//
real Complex:: Magnitude() const {

  return sqrt(Real*Real+Imag*Imag);
}

real Complex:: Power() const {

  return Real*Real+Imag*Imag;
}

real Complex:: PhaseB() const {
  if(Real!=.0 || Imag!=.0)
    return atan(Imag/Real);
  else
    return .0;
}

real Complex:: Phase() const {
  if(Real!=.0 || Imag!=.0)
    return atan2(Imag, Real);
  else
    return .0;
}

//
// IO operators
//
istream& operator>>(istream& s, Complex& comp_num){
  real val;

  s >> val;
  comp_num.SetReal(val);
  s >> val;
  comp_num.SetImag(val);

  return s;
}
ostream& operator<<(ostream& s, const Complex& comp_num) {
  s << "(" << comp_num.GetReal() << ",";
  s << comp_num.GetImag() << ")";

  return s;
}

