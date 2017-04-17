#include <CoreParticleUtil.h>
#include <CoreParticleUtilTemplate.h>
#include <DisjointSet.h>
#include <mex.h>
#include <map>

bool
surfaceParticle(CoreParticle* p, int ndim, bool bEightNeighbor)
{
	vector<vector<int>> nbh = NeighborhoodFactory::getInstance(ndim).neighbor4;
	if (bEightNeighbor)
	{
		nbh = NeighborhoodFactory::getInstance(ndim).neighbor8;
	}
	return p->neighbors.size() < nbh.size();
}

bool
medialParticle(CoreParticle* p, int ndim)
{
	return p->ascendents.size() >= 2 * (ndim - 1) || p->descendents.size() == 0;
}

bool
inflectionParticle(CoreParticle* p, int ndim)
{
	return surfaceParticle(p, ndim) && p->descendents.size() >= ndim;
}


CoreParticle*
coreParticleNdim(vector<int>& loc, int ndim)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	switch (ndim)
	{
	case 1:
		return factory.makeParticle(loc[0]);
	case 2:
		return factory.makeParticle(loc[0], loc[1]);
	case 3:
		return factory.makeParticle(loc[0], loc[1], loc[2]);
	case 4:
		return factory.makeParticle(loc[0], loc[1], loc[2], loc[3]);
	default:
		mexErrMsgTxt("particleNdim: unsupported number of dimensions. It has to be between 1 and 4.");
		return NULL;
	}
}

vector<int>
coreParticle2Index(CoreParticle*  p, int ndim)
{
	vector<int> idx(ndim);
	if (ndim == 1)
	{
		idx[0] = p->x;
	}
	else if (ndim == 2)
	{
		idx[0] = p->x;
		idx[1] = p->y;
	}
	else if (ndim == 3)
	{
		idx[0] = p->x;
		idx[1] = p->y;
		idx[2] = p->z;
	}
	else if (ndim == 4)
	{
		idx[0] = p->x;
		idx[1] = p->y;
		idx[2] = p->z;
		idx[3] = p->t;
	}
	else
	{
		mexErrMsgTxt("SetVoxel: unsupported number of dimensions. It has to be between 1 and 4.");
	}
	return idx;
}


vector<CoreParticle*>
generateParticleMap(vector<unsigned char>& L,
int ndim,
const int* dims)
{
	CoreParticleFactory& factory = CoreParticleFactory::getInstance();
	vector<CoreParticle*> mp(L.size(), NULL); //map to resolve uniqueness at each pixel
	for (int i = 0; i < L.size(); ++i)
	{
		if (L[i])
		{
			vector<int> sub = Ind2Sub(i, ndim, dims);
			CoreParticle* p = coreParticleNdim(sub, ndim);
			mp[i] = p;
		}
	}
	return mp;
}

vector<CoreParticle*>
setupParticleNeighbors(vector<CoreParticle*>& mp,
int ndim,
const int* dims,
bool bEightNeighbor)
{
	vector<CoreParticle*> particles;
	for (int i = 0; i < mp.size(); ++i)
	{
		if (mp[i] != NULL)
		{
			particles.push_back(mp[i]);
		}
	}
	NeighborhoodFactory& nfactory = NeighborhoodFactory::getInstance(ndim);
	vector<vector<int>> nbh = nfactory.neighbor4;
	if (bEightNeighbor)
	{
		nbh = nfactory.neighbor8;
	}

	for (int i = 0; i < particles.size(); ++i)
	{
		CoreParticle* p = particles[i];
		vector<int> sub = coreParticle2Index(p, ndim);
		for (int n = 0; n < nbh.size(); ++n)
		{
			if (NeighborCheck(sub.begin(), nbh[n].begin(), ndim, dims))
			{
				int idx = Sub2Ind(sub, nbh[n], ndim, dims);
				if (mp[idx] != NULL)
				{
					p->neighbors.push_back(mp[idx]);
				}
			}
		}
	}
	return particles;
}

/*
Cluster particles using disjoint set.
*/
vector<vector<CoreParticle*>>
clusterParticles(vector<CoreParticle*>& particles)
{
	set<CoreParticle*> S;
	vector<Node<CoreParticle*>*> nodes;
	map<CoreParticle*, int> imap;
	for (int i = 0; i < particles.size(); ++i)
	{
		Node<CoreParticle*>* n = makeset(particles[i]);
		nodes.push_back(n);
		imap[particles[i]] = i;
		S.insert(particles[i]);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		for (int j = 0; j < p->neighbors.size(); ++j)
		{
			CoreParticle* q = p->neighbors[j];
			if (S.find(q) != S.end())
			{
				merge(nodes[i], nodes[imap[q]]);
			}
		}
	}
	vector<Node<CoreParticle*>*> rep = clusters(nodes);
	vector<vector<CoreParticle*>> group(rep.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = distance(rep.begin(), find(rep.begin(), rep.end(), findset(nodes[i])));
		group[k].push_back(nodes[i]->key);
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return group;
}

CoreParticle*
representativeDescendent(CoreParticle* p)
{
	float dx = 0, dy = 0, dz = 0, dt = 0;
	int count = 0;
	for (int i = 0; i < p->neighbors.size(); ++i)
	{
		CoreParticle* q = p->neighbors[i];
		if (p->gen < q->gen)
		{
			dx += (q->x - p->x);
			dy += (q->y - p->y);
			dz += (q->z - p->z);
			dt += (q->t - p->t);
			count++;
		}
	}
	if (count == 0) return NULL;
	else
	{
		float mx = p->x + dx / count;
		float my = p->y + dy / count;
		float mz = p->z + dz / count;
		float mt = p->t + dt / count;
		CoreParticle* best = NULL;
		float mind = std::numeric_limits<float>::infinity();
		for (int i = 0; i < p->neighbors.size(); ++i)
		{
			CoreParticle* q = p->neighbors[i];
			if (p->gen < q->gen)
			{
				float d = (q->x - mx)*(q->x - mx) + (q->y - my)*(q->y - my) + (q->z - mz)*(q->z - mz) + (q->t - mt)*(q->t - mt);
				if (d < mind)
				{
					mind = d;
					best = q;
				}
			}
		}
		return best;
	}
}