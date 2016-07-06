#ifndef ___SZ_ISOTROPY_H___
#define ___SZ_ISOTROPY_H___

#include<vector>
#include<algorithm>
using namespace std;

#include "szMexUtility.h"

void
IsotropyMeasure(vector<float>& D, 
                vector<unsigned char>& L,
                vector<float>& vWall, 
                vector<float>& vVessel, 
                vector<float>& vNodule, 
                const int* dims);

#endif /* ___SZ_ISOTROPY_H___ */