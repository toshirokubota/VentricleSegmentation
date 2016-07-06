#include <mex.h>
#include <szRetrace.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>
#include <szParticle.h>

int
traceUp(vector<int>& M, 
		const vector<float>& D, 
		const CParticle& p, 
		const vector<CParticle>& vppeaks, 
		const int* dims)
{
	int x = p.m_X;
	int y = p.m_Y;
	int z = p.m_Z;
	int mval = GetData3(M, x, y, z, dims[0], dims[1], dims[2], (int)-1);
	if(mval >= 0)
	{
		return mval;
	}

	float dval1 = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	if(dval1 <= 0)
	{
		return -1;
	}

	//check if we reached a peak (or the volume boundary)
	vector<CParticle>::const_iterator pi = find(vppeaks.begin(), vppeaks.end(), p);
	int id;
	if(pi < vppeaks.end())
	{
		id = (int)(pi - vppeaks.begin());
		SetData3(M, x, y, z, dims[0], dims[1], dims[2], id);
		return id;
	}

	float dmax = 0;
	vector<CParticle> pn;
	for(int i=0; i<NumNeighbors; ++i)
	{
		int x2 = x+XOffset[i];
		int y2 = y+YOffset[i];
		int z2 = z+ZOffset[i];
		float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
		if(dval2 <= 0)
		{
			continue;
		}
		else if(dmax < dval2)
		{
			pn.clear();
			pn.push_back(CParticle(x2, y2, z2));
			dmax = dval2;
		}
		else if(dmax == dval2)
		{
			pn.push_back(CParticle(x2, y2, z2));
		}
	}
	if(dmax > dval1)
	{
		float dist = -1;
		int jd = -1;
		for(int j=0; j<pn.size(); ++j)
		{
			int id = traceUp(M, D, pn[j], vppeaks, dims);
			if(id >= 0)
			{
				int x2 = vppeaks[id].m_X;
				int y2 = vppeaks[id].m_Y;
				int z2 = vppeaks[id].m_Z;
				float d = (x2-x)*(x2-x) + (y2-y)*(y2-y) + (z2-z)*(z2-z);
				if(dist < 0 || d < dist)
				{
					dist = d;
					jd = id;
				}
			}
		}
		if(jd >= 0)
		{
			SetData3(M, x, y, z, dims[0], dims[1], dims[2], jd);
			return jd;
		}
	}
	else //the current position is on a local maximum
	{
		//this shouldn't happen, as all peaks are spotted at the beginning.
		return -1;
	}
	return -1;
}

vector<CParticle>
collectTraced(const vector<int>& M,
			  const vector<CParticle>& vppeaks,
			  const int* dims)
{
	vector<CParticle> vseeds;
	int i, j, k;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				int m = GetData3(M, k, j, i, dims[0], dims[1], dims[2], (int)-1);
				if(m>=0 && m<vppeaks.size())
				{
					if(find(vseeds.begin(), vseeds.end(), vppeaks[m]) == vseeds.end())
					{
						vseeds.push_back(vppeaks[m]);
					}
				}
			}
		}
	}
	return vseeds;
}
