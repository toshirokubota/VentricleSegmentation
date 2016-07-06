#include <mex.h>
#include <szTrace.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>

//#define STRICT_DOWNHILL

void
depthFirstTrace(vector<unsigned char>& B,
				const vector<float>& D,
				int x, int y, int z,
				const int* dims)
{
	if(GetData3(B, x, y, z, dims[0], dims[1], dims[2], (unsigned char)1))
	{
		return;
	}
	if(GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0) <= 0)
	{
		return;
	}

	SetData3(B, x, y, z, dims[0], dims[1], dims[2], (unsigned char)1);
	float dval = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	float mind = dval;
	vector<int> vx;
	vector<int> vy;
	vector<int> vz;
	for(int j=0; j<NumNeighbors; ++j)
	{
		int x2 = x + XOffset[j];
		int y2 = y + YOffset[j];
		int z2 = z + ZOffset[j];
		float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
		if(dval2 < mind)
		{
			mind = dval2;
			vx = vector<int>(1, x2);
			vy = vector<int>(1, y2);
			vz = vector<int>(1, z2);
		}
		else if(dval2 == mind)
		{
			vx.push_back(x2);
			vy.push_back(y2);
			vz.push_back(z2);
		}
	}
	for(int j = 0; j<vx.size(); ++j)
	{
		depthFirstTrace(B, D, vx[j], vy[j], vz[j], dims);
	}
}


void
SampledTrace(vector<unsigned char>& B,
			 const vector<float>& D,
			 const vector<CParticle>& vSeeds,
			 const int* dims)
{
	int i, j;
	for(i=0; i<vSeeds.size(); ++i)
	{
		int x = vSeeds[i].m_X;
		int y = vSeeds[i].m_Y;
		int z = vSeeds[i].m_Z;
		float dval = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
		for(j=0; j<NumNeighbors; ++j)
		{
			int x2 = x + XOffset[j];
			int y2 = y + YOffset[j];
			int z2 = z + ZOffset[j];
			float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
			if(dval2 < dval)
			{
				depthFirstTrace(B, D, x2, y2, z2, dims);
			}
		}
	}
}

void
incrementalTrace(vector<unsigned char>& B,
				 const vector<float>& D,
				 const vector<CParticle>& vSeeds,
				 const int* dims)
{
	vector<CParticle> vp = vSeeds;
	int initDvalSq = 0;
	int i, j, k;
	for(int i=0; i<vp.size(); ++i)
	{
		int x1=vp[i].m_X;
		int y1=vp[i].m_Y;
		int z1=vp[i].m_Z;
		SetData3(B, x1, y1, z1, dims[0], dims[1], dims[2], (unsigned char)ForegroundColor);
		float dval = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
		int ivalSq = (int)(dval*dval);
		if(ivalSq > initDvalSq)
		{
			initDvalSq = ivalSq;
		}
	}
	printf("incrementalTrace: init dvalue = %d\n", initDvalSq);
	for(int idval = initDvalSq; idval > 0; idval--)
	{
		vector<CParticle> vp;
		for(i=0; i<dims[2]; ++i) 
		{ 
			for(j=0; j<dims[1]; ++j) 
			{ 
				for(k=0; k<dims[0]; ++k) 
				{ 
					if(GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)1) == 0)
					{
						float dval = GetData3(D, k, j, i, dims[0], dims[1], dims[2], (float)0);
						int idvalSq = (int)(dval*dval);
						if(idvalSq == idval)
						{
							if(onBoundary3(B, k, j, i, dims))
							{
								vp.push_back(CParticle(k, j, i));
							}
						}
					}
				}
			}
		}
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			SetData3(B, x1, y1, z1, dims[0], dims[1], dims[2], (unsigned char)ForegroundColor);
		}
	}
}

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
		//printf("Free-growth: #vp = %d, #vp2 = %d\n", vp.size(), vp2.size());
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
					if(d < dthres)
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
#ifdef STRICT_DOWNHILL
						if(dval1 > dval2)
#else
						if(dval1 >= dval2)
#endif
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

/*void
ConstrainedGrowth(vector<unsigned char>& B,
				  const vector<float>& D,
				  vector<unsigned char>& P,
				  float dthres,
				  int x, int y, int z,
				  const int* dims)
{
	const int PosLabel = 1; //label for a nodule
	const int NegLabel = 2; //label for a non-nodule

	int i, j, k;
	vector<CParticle> vpP;
	vector<CParticle> vpN;

	float xf = x, yf = y, zf = z;
	vpP.push_back(CParticle(x, y, z, PosLabel));
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					float d = sqrt((float)(xf-k)*(xf-k)+(yf-j)*(yf-j)+(zf-i)*(zf-i));
					if(d < dthres)
					{
						vpP.push_back(CParticle(k, j, i, PosLabel));
					}
					else
					{
						vpN.push_back(CParticle(k, j, i, NegLabel));
					}
				}
			}
		}
	}
	vector<unsigned char> Bn(B.size(), 0); 

	FreeGrowth(B, D, vpP, dims);		
	FreeGrowth(Bn, D, vpN, dims);		

	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(Bn, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0) &&
					GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					//SetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)128);
					float mindsp = -1.0, mindsn = -1.0;
					for(int m=0; m<vpP.size(); ++m)
					{
						int x1 = vpP[m].m_X;
						int y1 = vpP[m].m_Y;
						int z1 = vpP[m].m_Z;
						float dsq = (x1-k)*(x1-k) + (y1-j)*(y1-j) + (z1-i)*(z1-i);
						if(mindsp < 0 || dsq < mindsp)
						{
							mindsp = dsq;
						}
					}
					for(int m=0; m<vpN.size(); ++m)
					{
						int x1 = vpN[m].m_X;
						int y1 = vpN[m].m_Y;
						int z1 = vpN[m].m_Z;
						float dsq = (x1-k)*(x1-k) + (y1-j)*(y1-j) + (z1-i)*(z1-i);
						if(mindsn < 0 || dsq < mindsn)
						{
							mindsn = dsq;
						}
					}
					if(mindsp > mindsn)
					{
						SetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0);
					}
				}
			}
		}
	}
}
*/
