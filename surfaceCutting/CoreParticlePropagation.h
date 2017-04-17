#ifndef __CORE_PARTICLE_PROPAGATION_H__
#define __CORE_PARTICLE_PROPAGATION_H__

#include <CoreParticle.h>
#include <Graph.h>
#include <GraphFactory.h>

/*
from each surface particle, move towrad the center and establish ascendent/descendent relation.
*/
vector<CoreParticle*>
propagateParticles(vector<CoreParticle*>& particles,
					vector<CoreParticle*>& mp,
					vector<int>& S,
					int ndim, const int* dims);


vector<Edge<CoreParticle*>*>
makeGraphStructure(vector<CoreParticle*>& mp, 
					vector<int>& S, 
					int ndim, const int* dims);

#endif /* __CORE_PARTICLE_PROPAGATION_H__ */