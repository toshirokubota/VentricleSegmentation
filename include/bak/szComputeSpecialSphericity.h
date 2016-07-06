#ifndef _SZ_COMPUTE_SPECIAL_SPHERICITY_H_
#define _SZ_COMPUTE_SPECIAL_SPHERICITY_H_

#include <vector>
using namespace std;

void
ComputeSphericity(vector<float>& S, 
                  vector<float>& V, 
                  vector<unsigned char>& L, 
                  float idealV,
                  int w, 
                  const int* dims);

#endif /* _SZ_COMPUTE_SPECIAL_SPHERICITY_H_ */