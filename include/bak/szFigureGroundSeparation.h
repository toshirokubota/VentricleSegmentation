#ifndef _SZ_FIGURE_GROUND_SEPARATION_H_
#define _SZ_FIGURE_GROUND_SEPARATION_H_

#include<vector>
using namespace std;

void
ReactionDiffusion(vector<float>& A, 
                  int xD, 
                  int yD, 
                  int zD, 
                  float lambda,
				  float dfspeed[6],
				  float reacWeight);


bool
FigureGroundSeparation(const vector<unsigned short>& V,  //input sub-volume
                       vector<unsigned char>& L,  //output segmentation
					   vector<float>& Q, //output probablity image
                       int mode,                  //0: solid, 1: non-solid
                       int numiter,
                       float lambda,
                       const int* dimsA);
  

#endif /* _SZ_FIGURE_GROUND_SEPARATION_H_ */