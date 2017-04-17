#ifndef __CORE_PARTICLE_UTIL_H__
#define __CORE_PARTICLE_UTIL_H__

#include <CoreParticle.h>
#include <NeighborhoodND.h>

CoreParticle*
coreParticleNdim(vector<int>& loc, int ndim);

vector<int>
coreParticle2Index(CoreParticle*  p, int ndim);

bool
surfaceParticle(CoreParticle* p, int ndim, bool bEightNeighbor = false);

bool
medialParticle(CoreParticle* p, int ndim);

bool
inflectionParticle(CoreParticle* p, int ndim);

vector<CoreParticle*>
generateParticleMap(vector<unsigned char>& L, int ndim, const int* dims);

vector<CoreParticle*>
setupParticleNeighbors(vector<CoreParticle*>& mp, int ndim, const int* dims, bool bEightNeighbor=false);

vector<vector<CoreParticle*>>
clusterParticles(vector<CoreParticle*>& particles);

CoreParticle*
representativeDescendent(CoreParticle* p);

#endif /* __CORE_PARTICLE_UTIL_H__ */