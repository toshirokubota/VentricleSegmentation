#ifndef _SZ_IS_FUNCTIONS_H_
#define _SZ_IS_FUNCTIONS_H_

#include <vector>
using namespace std;

bool
IsDegenerate(const vector<unsigned char>& L);

bool
IsSolitary(const vector<unsigned char>& L, 
           int max_size, 
           const int* dims);

bool
IsSolitary(const vector<unsigned char>& L, 
           int x, int y, int z,
           int max_size, 
           const int* dims);

bool
IsWallAttached(const vector<unsigned char>& L,
               const int* dims);

#endif /* _SZ_IS_FUNCTIONS_H_ */