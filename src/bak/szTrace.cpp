#include <mex.h>
#include <szTrace.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szMyNeighborOp.h>

#define STRICT_DOWNHILL

void
FreeGrowth(vector<unsigned char>& B,
		   const vector<float>& D,
		   const vector<CParticle>& vSeeds,
		   const int* dims)
{
	vector<CParticle> vp = vSeeds;
	int i, j;
	while(!vp.empty())
	{
		vector<CParticle> vp2;
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			SetData3(B, x1, y1, z1, dims[0], dims[1], dims[2], (unsigned char)ForegroundColor);
			float dval1 = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
			for(j=0; j<NumNeighbors; ++j)
			{
				int x2 = x1+XOffset[j];
				int y2 = y1+YOffset[j];
				int z2 = z1+ZOffset[j];
				if(!GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)1))
				{
					float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
					if(dval2>0)
					{
#ifdef STRICT_DOWNHILL
						if(dval1 > dval2)
#else
						if(dval1 >= dval2)
#endif
						{
							CParticle pn(x2, y2, z2, 1);
							if(find(vp2.begin(), vp2.end(), pn) == vp2.end())
							{
								vp2.push_back(pn);
							}
						}
					}
				}
			}
		}
		printf("Free-growth: #vp = %d, #vp2 = %d\n", vp.size(), vp2.size());
		vp = vp2;
	}
}

void
ConstrainedGrowth(vector<unsigned char>& B,
				  const vector<float>& D,
				  vector<unsigned char>& P,
				  int x, int y, int z,
				  const int* dims)
{
	const int PosLabel = 1; //label for a nodule
	const int NegLabel = 2; //label for a non-nodule

	int i, j, k;
	vector<CParticle> vp;
	vector<int> C(B.size(), 0); 
	/*int numClusters = ConnectedComponentAnalysisBigger(C, P, NeighborhoodFour, (unsigned char)0, 3, dims);
	int label = GetData3(C, x, y, z, dims[0], dims[1], dims[2], 0);
	printf("ConstrainedGrowth: label = %d\n", label);
	float xf = 0, yf = 0, zf = 0;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(C, k, j, i, dims[0], dims[1], dims[2], -1) == label)
				{
					vp.push_back(CParticle(k, j, i, PosLabel));
					SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)PosLabel);
					xf += k;
					yf += j;
					zf += i;
				}
			}
		}
	}
	xf /= vp.size();
	yf /= vp.size();
	zf /= vp.size();*/

	float xf = x, yf = y, zf = z;
	vp.push_back(CParticle(x, y, z, PosLabel));
	float dthres = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	printf("ConstrainedGrowth: (x, y, z) = (%d, %d, %d), dthres = %f\n", 
		x+1, y+1, z+1, dthres);
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					float d = sqrt((float)(xf-k)*(xf-k)+(yf-j)*(yf-j)+(zf-i)*(zf-i));
					if(d <= dthres)
					{
						vp.push_back(CParticle(k, j, i, PosLabel));
						SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)PosLabel);
						//printf("ConstrainedGrowth: Setting positive\n");
					}
					else
					{
						vp.push_back(CParticle(k, j, i, NegLabel));
						SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)NegLabel);
						//printf("ConstrainedGrowth: Setting negative\n");
					}
				}
			}
		}
	}

	int time_stamp = 0;
	while(!vp.empty())
	{
		time_stamp ++;
		vector<CParticle> vp2;
		vector<CParticle> vrp;
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			SetData3(B, x1, y1, z1, dims[0], dims[1], dims[2], (unsigned char)vp[i].m_Life);
			float dval1 = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
			for(j=0; j<NumNeighbors; ++j)
			{
				int x2 = x1+XOffset[j];
				int y2 = y1+YOffset[j];
				int z2 = z1+ZOffset[j];
				if(!GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)1))
				{
					float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
					if(dval2>0)
					{
						if(dval1 >= dval2)
						{
							CParticle pn(x2, y2, z2, vp[i].m_Life);
							if(find(vp2.begin(), vp2.end(), pn) == vp2.end())
							{
								vp2.push_back(pn);
							}
						}
					}
				}
			}
		}
		printf("Constrained-growth: #vp = %d, #vp2 = %d\n", vp.size(), vp2.size());
		vp = vp2;
	}
	for(int i=0; i<B.size(); ++i)
	{
		if(B[i] == PosLabel)
			B[i] = ForegroundColor;
		else
			B[i] = 0;
	}
}

