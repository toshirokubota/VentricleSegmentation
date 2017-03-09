#ifndef _SZ_REACTION_DIFFUSION_H_
#define _SZ_REACTION_DIFFUSION_H_

#include<vector>
using namespace std;

void
ReactionDiffusion3D(
					vector<float>& A,
					int xD,
					int yD,
					int zD,
					vector<float>& C,
					float dfspeed[6],
					float rate);

void
ReactionDiffusion4D(
					vector<float>& A,
					int xD,
					int yD,
					int zD,
					int tD,
					vector<float>& C,
					float dfspeed[8],
					float rate);


#endif /* _SZ_REACTION_DIFFUSION_H_ */