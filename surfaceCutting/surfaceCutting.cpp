#include <iostream>
using namespace std;
#include <stdio.h>

#ifdef MEX_DLL
#include <mex.h>
#include "mexFileIO.h"
#endif

#include <vector>
#include <algorithm>
#include <set>
#include <map>
using namespace std;
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <szDefaultParam.h>
#include <szMiscOperations.h>
#include <DisjointSet.h>
#include <szParticle4D.h>
#include <szConnectedComponent.h>

#include <CoreParticle.h>
#include <NeighborhoodND.h>
#include <CoreParticleUtil.h>

struct CuttingPlane {
	CParticle4D orgin;
	vector<CParticle4D> adjacents;
};

float
dotProduct(float x1, float y1, float z1, float t1, float x2, float y2, float z2, float t2)
{
	return x1*x2 + y1*y2 + z1*z2 + t1*t2;
}

float
length(float x, float y, float z, float t)
{
	return sqrt(x*x + y*y + z*z + t*t);
}

vector<int> cutSurfaceByPlane(vector<CoreParticle*>& surface, CuttingPlane plane, int ndim)
{
	CParticle4D orig = plane.orgin;
	vector<int> label(surface.size(), 0);
	vector<vector<float>> dirs;
	for (int i = 0; i < plane.adjacents.size(); ++i)
	{
		vector<float> dir(4, 0.0f);
		dir[0] = plane.adjacents[i].m_X - orig.m_X;
		dir[1] = plane.adjacents[i].m_Y - orig.m_Y;
		dir[2] = plane.adjacents[i].m_Z - orig.m_Z;
		dir[3] = plane.adjacents[i].m_T - orig.m_T;
		dirs.push_back(dir);
	}
	vector<float> perp(4, 0.0f);
	if (ndim == 2)
	{
		perp[0] = dirs[0][1];
		perp[1] = -dirs[0][0];
		if (perp[0] * dirs[1][0] + perp[1] * dirs[1][1] > 0)
		{
			perp[0] = -perp[0];
			perp[1] = -perp[1];
		}
	}
	else if (ndim == 3)
	{
		vector<float> a = dirs[0];
		vector<float> b = dirs[1];
		perp[0] = a[1] * b[2] - a[2] * b[1];
		perp[1] = a[2] * b[0] - a[0] * b[2];
		perp[2] = a[0] * b[1] - a[1] * b[0];
		if (perp[0] * dirs[2][0] + perp[1] * dirs[2][1] + perp[2] * dirs[2][2]> 0)
		{
			perp[0] = -perp[0];
			perp[1] = -perp[1];
			perp[2] = -perp[2];
		}
	}
	else
	{
		mexErrMsgTxt("surfaceCutting: Only 2-dim and 3-dim are supported at this time.");
	}

	//separate points by the cutting plane
	vector<Node<CoreParticle*>*> nodes;
	map<CoreParticle*, Node<CoreParticle*>*> nmap;
	for (int i = 0; i < surface.size(); ++i)
	{
		CoreParticle* p = surface[i];
		float dp = dotProduct(perp[0], perp[1], perp[2], perp[3], p->x - orig.m_X, p->y - orig.m_Y, p->z - orig.m_Z, p->t - orig.m_T);
		p->label = dp > 0 ? 1 : (dp < 0 ? -1 : 0);
		Node<CoreParticle*>* n = makeset(p);
		nodes.push_back(n);
		nmap[p] = n;
	}
	//form clusters
	for (int i = 0; i < nodes.size(); ++i)
	{
		CoreParticle* p = nodes[i]->key;
		for (int j = 0; j < p->neighbors.size(); ++j)
		{
			CoreParticle* q = p->neighbors[j];
			if (p->label == q->label && nmap.find(q) != nmap.end())
			{
				merge(nodes[i], nmap[q]);
			}
		}
	}
	vector<Node<CoreParticle*>*> rep = clusters(nodes);
	map<Node<CoreParticle*>*, int> imap;
	for (int i = 0; i < rep.size(); ++i)
	{
		imap[rep[i]] = i;
	}
	vector<vector<Node<CoreParticle*>*>> group(rep.size());
	for (int i = 0; i < nodes.size(); ++i)
	{
		int k = imap[findset(nodes[i])];
		group[k].push_back(nodes[i]);
	}
	//find a cluster with the origin
	int ocluster = -1;
	float mind = std::numeric_limits<float>::infinity();
	for (int i = 0; i < rep.size(); ++i)
	{
		if (rep[i]->key->label == 0)
		{
			for (int j = 0; j < group[i].size(); ++j)
			{
				CoreParticle* p = group[i][j]->key;
				float len = length(p->x - orig.m_X, p->y - orig.m_Y, p->z - orig.m_Z, p->t - orig.m_T);
				if (len < mind)
				{
					mind = len;
					ocluster = i;
				}
			}
		}
	}

	/*
	Among clusters with non-zero label, find the one based on the following criteria:
	1. it has to be adjacent to the cluster containing the origin.
	2. most particles project to the positive side of the project plane.
	*/
	//select the one that is closest to one of adjacent point in the plane
	vector<int> cand;
	for (int i = 0; i < rep.size(); ++i)
	{
		if (rep[i]->key->label != 1) continue; //only the label of 1 is relevant
		for (int j = 0; j<group[i].size(); ++j)
		{
			CoreParticle* p = group[i][j]->key;
			bool done = false;
			for (int k = 0; k < group[ocluster].size(); ++k)
			{
				CoreParticle* q = group[ocluster][k]->key;
				float d = length(p->x - q->x, p->y - q->y, p->z - q->z, p->t - q->t);
				if (d <= 1.0f)
				{
					done = true;
					break;
				}
			}
			if (done)
			{
				cand.push_back(i);
				break;
			}
		}
	}
	//now select the cluster
	Node<CoreParticle*>* chosen = NULL;
	if (cand.size() == 1)
	{
		chosen = rep[cand[0]];
	}
	else
	{
		float maxscore = 0;
		int best = -1;
		for (int c = 0; c < cand.size(); ++c)
		{
			int count = 0; 
			int k = cand[c];
			for (int i = 0; i < group[k].size(); ++i)
			{
				CoreParticle* p = group[k][i]->key;
				for (int j = 0; j < ndim - 1; ++j)
				{
					float dp = dotProduct(p->x - orig.m_X, p->y - orig.m_Y, p->z - orig.m_Z, p->t - orig.m_T,
						dirs[j][0], dirs[j][1], dirs[j][2], dirs[j][3]);
					if (dp > 0) count++;
				}
			}
			float score = (float)count / group[k].size();
			if (score > maxscore)
			{
				maxscore = score;
				best = k;
			}
		}
		chosen = rep[best];
	}
	for (int i = 0; i < nodes.size(); ++i)
	{
		if (findset(nodes[i]) != chosen)
		{
			nodes[i]->key->label = 0;
		}
	}
	//also include the origin
	for (int i = 0; i < group[ocluster].size(); ++i)
	{
		CoreParticle* p = group[ocluster][i]->key;
		float sum = 0;
		for (int j = 0; j < ndim - 1; ++j)
		{
			float dp = dotProduct(p->x - orig.m_X, p->y - orig.m_Y, p->z - orig.m_Z, p->t - orig.m_T,
				dirs[j][0], dirs[j][1], dirs[j][2], dirs[j][3]);
			sum += dp;
		}
		//but only if it is on the positive side of the projection axis.
		if (sum > 0)
		{
			group[ocluster][i]->key->label = 1;
		}
	}
	//set the label output
	for (int i = 0; i < surface.size(); ++i)
	{
		if (surface[i]->label)
		{
			label[i] = 1;
		}
	}

	for (int i = 0; i < nodes.size(); ++i)
	{
		delete nodes[i];
	}
	return label;
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 0)
	{
		mexErrMsgTxt("Usage: [S] = particleGrowth(L, [thres mincore])");
		return;
	}

	//load figure-ground segmentation
	int ndimL;
	const int* dimsL;
	mxClassID classL;
	vector<unsigned char> L;
	LoadData(L, prhs[0], classL, ndimL, &dimsL);

	//Load cutting planes
	vector<CuttingPlane> planes;
	{
		int ndimD;
		const int* dimsD;
		mxClassID classD;
		vector<int> D;
		LoadData(D, prhs[1], classD, ndimD, &dimsD);
		for (int i = 0; i < dimsD[0]; ++i) {
			CuttingPlane plane;
			plane.orgin.m_X = GetData2(D, i, 0, dimsD[0], dimsD[1], 0);
			plane.orgin.m_Y = GetData2(D, i, 1, dimsD[0], dimsD[1], 0);
			plane.orgin.m_Z = GetData2(D, i, 2, dimsD[0], dimsD[1], 0);
			plane.orgin.m_T = GetData2(D, i, 3, dimsD[0], dimsD[1], 0);
			for (int j = 4; j < dimsD[1]; j+=4)
			{
				CParticle4D pc;
				pc.m_X = GetData2(D, i, j, dimsD[0], dimsD[1], 0);
				pc.m_Y = GetData2(D, i, j + 1, dimsD[0], dimsD[1], 0);
				pc.m_Z = GetData2(D, i, j + 2, dimsD[0], dimsD[1], 0);
				pc.m_T = GetData2(D, i, j + 3, dimsD[0], dimsD[1], 0);
				plane.adjacents.push_back(pc);
			}
			planes.push_back(plane);
		}
	}
	bool bEightNeighbor = true; //8-neighbor for cutting
	/*if (nrhs >= 3)
	{
		mxClassID classMode;
		int value;
		ReadScalar(value, prhs[2], classMode);
		bEightNeighbor = value > 0 ? true : false;
	}*/


	int nvoxels = numberOfElements(ndimL, dimsL);

	vector<CoreParticle*> mp = generateParticleMap(L, ndimL, dimsL);
	vector<CoreParticle*> particles = setupParticleNeighbors(mp, ndimL, dimsL, bEightNeighbor);
	vector<CoreParticle*> surface;
	for (int i = 0; i < particles.size(); ++i)
	{
		if (surfaceParticle(particles[i], ndimL, bEightNeighbor))
		{
			surface.push_back(particles[i]);
		}
	}
	vector<vector<int>> labels;
	for (int c = 0; c < planes.size(); ++c)
	{
		vector<int> lb = cutSurfaceByPlane(surface, planes[c], ndimL);
		labels.push_back(lb);
	}

	if (nlhs >= 1)
	{
		int dims[] = { surface.size(), 4};
		vector<int> F(dims[0] * dims[1]);
		for (int i = 0; i < surface.size(); ++i)
		{
			CoreParticle* p = surface[i];
			SetData2(F, i, 0, dims[0], dims[1], p->x);
			SetData2(F, i, 1, dims[0], dims[1], p->y);
			SetData2(F, i, 2, dims[0], dims[1], p->z);
			SetData2(F, i, 3, dims[0], dims[1], p->t);
		}
		plhs[0] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	if (nlhs >= 2)
	{
		int dims[] = { surface.size(), labels.size() };
		vector<int> F(dims[0] * dims[1]);
		int i = 0;
		for (int j = 0; j < labels.size(); ++j)
		{
			for (int i = 0; i < surface.size(); ++i)
			{
				SetData2(F, i, j, dims[0], dims[1], labels[j][i]);
			}
		}
		plhs[1] = StoreData(F, mxINT32_CLASS, 2, dims);
	}
	CoreParticleFactory::getInstance().clean();
	mexUnlock();
}

