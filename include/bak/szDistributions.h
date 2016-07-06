#ifndef ___SZ_DISTRIBUTIONS_H___
#define ___SZ_DISTRIBUTIONS_H___

#include <vector>
using namespace std;

vector<float>
computePDF(const vector<int>& vhist);

vector<float>
computeCDF(const vector<int>& vhist);

vector<float>
computeCDF(const vector<float>& vpdf);

float 
computeKSD(const vector<float>& ca,
		   const vector<float>& cb);

float 
computeJSD(const vector<float>& pa,
		   const vector<float>& pb);

#endif /* ___SZ_DISTRIBUTIONS_H___ */