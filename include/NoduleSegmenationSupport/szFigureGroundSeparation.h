#ifndef _SZ_FIGURE_GROUND_SEPARATION_H_
#define _SZ_FIGURE_GROUND_SEPARATION_H_

#include<vector>
using namespace std;

void
ReactionDiffusionInitialize(const vector<float>& A,  //input sub-volume
					   vector<float>& Q, //output probability image
                       int mode,                  //0: solid, 1: non-solid
                       const int* dimsA);

void
ReactionDiffusion(vector<float>& A, 
                  int xD, 
                  int yD, 
                  int zD, 
                  float lambda,
				  float dfspeed[6],
				  float reacWeight);

void
ReactionDiffusionNew(vector<float>& A, 
					 int xD, 
					 int yD, 
					 int zD,
					 vector<float>& C,
					 float dfspeed[6],
					 float rate);

void
ReactionDiffusionJacobi(vector<float>& A, 
						int xD, 
						int yD, 
						int zD,
						vector<float>& C,
						float dfspeed[6],
						float rate);

bool
FigureGroundSeparation(const vector<float>& V,  //input sub-volume
                       vector<unsigned char>& L,  //output segmentation
					   vector<float>& Q, //output probablity image
                       int mode,                  //0: solid, 1: non-solid
                       int numiter,
                       float lambda,
                       const int* dimsA);
  
bool
FigureGroundSeparationMasked(const vector<float>& A,  //input sub-volume
							 vector<unsigned char>& L,  //output segmentation
							 vector<float>& Q, //output probability image
							 vector<char>& M,  //input masked image (1:foreground, -1: background, 0: uncertain)
							 int mode,                  //0: solid, 1: non-solid
							 int numiter,
							 float lambda,
							 const int* dimsA);

#endif /* _SZ_FIGURE_GROUND_SEPARATION_H_ */