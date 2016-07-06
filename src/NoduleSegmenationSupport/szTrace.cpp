#include <mex.h>
#include <szTrace.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szMyNeighborOp.h>
#include <szMiscOperations.h>

#define STRICT_DOWNHILL

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
		printf("incrementalTrace: dval = %f at (%d, %d, %d)\n", dval, x1+1, y1+1, z1+1);
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
					if(onSurface3(B, k, j, i, dims))
					{
						vp.push_back(CParticle(k, j, i));
					}
				}
			}
		}
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			for(j=0; j<NumNeighbors; ++j)
			{
				int x2 = x1+XOffset[j];
				int y2 = y1+YOffset[j];
				int z2 = z1+ZOffset[j];
				float dval1 = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
				if(GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)1) == 0)
				{
					float dval = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
					int idvalSq = (int)(dval*dval);
					//if(idvalSq >= idval && dval <= dval1)
					if(idvalSq == idval)
					{
						SetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)ForegroundColor);
						//SetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)idval);
					}
				}
			}
		}
	}
}

int
FreeGrowth(vector<int>& B,
		   const vector<float>& D,
		   const vector<CParticle>& vSeeds,
		   bool bStrict,
		   const int* dims)
{
	vector<CParticle> vp = vSeeds;
	int i, j;
	int cnt = 1;
	while(!vp.empty())
	{
		cnt++;
		vector<CParticle> vp2;
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			SetData3(B, x1, y1, z1, dims[0], dims[1], dims[2], (int)cnt);
		}
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			float dval1 = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
			for(j=0; j<NumNeighbors; ++j)
			{
				int x2 = x1+XOffset[j];
				int y2 = y1+YOffset[j];
				int z2 = z1+ZOffset[j];
				if(!GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (int)1))
				{
					float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
					//printf("%d (%d %d %d) %f -> (%d %d %d) %f\n", cnt, x1+1, y1+1, z1+1, dval1, x2+1, y2+1, z2+1, dval2);
					if(dval2>0)
					{
						if((bStrict && dval1 > dval2) || (!bStrict && dval1 >= dval2))
						{
							CParticle pn(x2, y2, z2, 1);
							if(find(vp2.begin(), vp2.end(), pn) == vp2.end())
							{
								//printf("pushing (%d %d %d)\n", x2+1, y2+1, z2+1);
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
	return cnt;
}


int
FreeGrowth6(vector<int>& B,
		   const vector<float>& D,
		   const vector<CParticle>& vSeeds,
		   bool bStrict,
		   const int* dims)
{
	vector<CParticle> vp = vSeeds;
	int i, j;
	int cnt = 1;
	while(!vp.empty())
	{
		cnt++;
		vector<CParticle> vp2;
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			SetData3(B, x1, y1, z1, dims[0], dims[1], dims[2], (int)cnt);
		}
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			float dval1 = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
			for(j=0; j<NumNeighbors6; ++j)
			{
				int x2 = x1+XOffset6[j];
				int y2 = y1+YOffset6[j];
				int z2 = z1+ZOffset6[j];
				if(!GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (int)1))
				{
					float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (float)0);
					//printf("%d (%d %d %d) %f -> (%d %d %d) %f\n", cnt, x1+1, y1+1, z1+1, dval1, x2+1, y2+1, z2+1, dval2);
					if(dval2>0)
					{
						if((bStrict && dval1 > dval2) || (!bStrict && dval1 >= dval2))
						{
							CParticle pn(x2, y2, z2, 1);
							if(find(vp2.begin(), vp2.end(), pn) == vp2.end())
							{
								//printf("pushing (%d %d %d)\n", x2+1, y2+1, z2+1);
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
	return cnt;
}

/*
grows only to downhill directions and the gradient vector and growth direction (current position - starting position)
forms an acute angle.
*/
void
FreeGrowthAngle0(vector<unsigned char>& B,
				const vector<float>& D,
				const CParticle& seed,
				const int* dims)
{
	vector<CParticle> vp(1, seed);
	int x0 = seed.m_X;
	int y0 = seed.m_Y;
	int z0 = seed.m_Z;
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
							float dvalW = GetData3(D, x2-1, y2, z2, dims[0], dims[1], dims[2], (float)dval2);
							float dvalE = GetData3(D, x2+1, y2, z2, dims[0], dims[1], dims[2], (float)dval2);
							float dvalN = GetData3(D, x2, y2-1, z2, dims[0], dims[1], dims[2], (float)dval2);
							float dvalS = GetData3(D, x2, y2+1, z2, dims[0], dims[1], dims[2], (float)dval2);
							float dvalB = GetData3(D, x2, y2, z2-1, dims[0], dims[1], dims[2], (float)dval2);
							float dvalF = GetData3(D, x2, y2, z2+1, dims[0], dims[1], dims[2], (float)dval2);
							float dx = dvalW-dvalE;
							float dy = dvalN-dvalS;
							float dz = dvalB-dvalF;
							float dmag = sqrt(dx*dx + dy*dy + dz*dz);
							float px = x2-x0;
							float py = y2-y0;
							float pz = z2-z0;
							float pmag = sqrt(px*px + py*py + pz*pz);
							float dp = (px*dx + py*dy + pz*dz)/pmag*dmag;
							if(dp > 3.141529/4.0)
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
		}
		//printf("Free-growth: #vp = %d, #vp2 = %d\n", vp.size(), vp2.size());
		vp = vp2;
	}
}

/*
grows only to downhill directions and the gradient vector and growth direction (current position - starting position)
forms an acute angle.
*/
void
FreeGrowthAngle(vector<unsigned char>& B,
				const vector<float>& D,
				const vector<CParticle>& vseed,
				const int* dims)
{
	for(int i=0; i<vseed.size(); ++i)
	{
		vector<unsigned char> S(B.size(), 0);
		FreeGrowthAngle0(S, D, vseed[i], dims);
		for(int j=0; j<S.size(); ++j)
		{
			if(S[j])
			{
				B[j] = ForegroundColor;
			}
		}
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
						//vp.push_back(CParticle(k, j, i, PosLabel));
						vp.insert(vp.begin(), CParticle(k, j, i, PosLabel));
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

/*
Curve out the free-growth segmentation by those grown from negative seeds.
*/
void
NegativeMigrateGrowth(vector<unsigned char>& B,
					  const vector<unsigned char>& S,
					  const vector<float>& D,
					  const vector<unsigned char>& P,
					  const vector<CParticle>& vseeds,
					  const int* dims)
{
	const int PosLabel = 1; //label for a nodule
	const int NegLabel = 2; //label for a non-nodule

	int i, j, k;
	vector<CParticle> vp;
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					if(GetData3(S, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
					{
						bool bFound = false;
						for(int m=0; m<vseeds.size(); ++m)
						{
							float xf = vseeds[m].m_X;
							float yf = vseeds[m].m_Y;
							float zf = vseeds[m].m_Z;
							float dthres = GetData3(D, (int)xf, (int)yf, (int)zf, dims[0], dims[1], dims[2], (float)0);
							float d = sqrt((float)(xf-k)*(xf-k)+(yf-j)*(yf-j)+(zf-i)*(zf-i));
							if(d < dthres)
							{
								bFound = true;
								break;
							}
						}
						if(bFound)
						{
							//vp.push_back(CParticle(k, j, i, PosLabel));
							//vp.insert(vp.begin(), CParticle(k, j, i, PosLabel));
							//SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)PosLabel);
							printf("NegativeMigrateGrowth: Setting positive @(%d, %d, %d)\n", k, j, i);
						}
						else
						{
							vp.push_back(CParticle(k, j, i, NegLabel));
							//SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)NegLabel);
							printf("NegativeMigrateGrowth: Setting negative @(%d, %d, %d)\n", k, j, i);
						}
					}
				}
			}
		}
	}

	while(!vp.empty())
	{
		vector<CParticle> vp2;
		for(i=0; i<vp.size(); ++i)
		{
			int x1=vp[i].m_X;
			int y1=vp[i].m_Y;
			int z1=vp[i].m_Z;
			SetData3(B, x1, y1, z1, dims[0], dims[1], dims[2], (unsigned char)1);
			float dval1 = GetData3(D, x1, y1, z1, dims[0], dims[1], dims[2], (float)0);
			for(j=0; j<NumNeighbors; ++j)
			{
				int x2 = x1+XOffset[j];
				int y2 = y1+YOffset[j];
				int z2 = z1+ZOffset[j];
				if(!GetData3(B, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)1))
				{
					if(GetData3(S, x2, y2, z2, dims[0], dims[1], dims[2], (unsigned char)0))
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
		}
		printf("NegativeMigrateGrowth: #vp = %d, #vp2 = %d\n", vp.size(), vp2.size());
		vp = vp2;
	}
	for(int i=0; i<S.size(); ++i)
	{
		if(S[i] && B[i] == 0)
			B[i] = ForegroundColor;
		else
			B[i] = 0;
	}
}

void
NewRegionPartition(vector<unsigned char>& B,
				   const vector<float>& D,
				   vector<unsigned char>& P,
				   int x, int y, int z,
				   const int* dims)
{
	int i, j, k;
	vector<CParticle> vp;
	vector<CParticle> vn;
	//vector<int> C(B.size(), 0); 

	float xf = x, yf = y, zf = z;
	float dthres = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					float d = sqrt((float)(xf-k)*(xf-k)+(yf-j)*(yf-j)+(zf-i)*(zf-i));
					float dval = GetData3(D, k, j, i, dims[0], dims[1], dims[2], (float)0);
					if(d < dthres)
					{
						vp.push_back(CParticle(k, j, i, dval));
						SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)3);
					}
					else
					{
						vn.push_back(CParticle(k, j, i, dval));
						SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)2);
					}
				}
			}
		}
	}

	vector<CParticle> vc = vp;

	int time_stamp = 0;
	while(!vc.empty())
	{
		time_stamp ++;
		vector<CParticle> vc2;
		for(i=0; i<vc.size(); ++i)
		{
			int x1=vc[i].m_X;
			int y1=vc[i].m_Y;
			int z1=vc[i].m_Z;
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
							CParticle pn(x2, y2, z2);
							if(find(vc2.begin(), vc2.end(), pn) == vc2.end())
							{
								float dminP = 1.0e10; //a big number
								for(int i2=0; i2<vp.size(); ++i2)
								{
									int x3 = vp[i2].m_X;
									int y3 = vp[i2].m_Y;
									int z3 = vp[i2].m_Z;
									float dval3 = vp[i2].m_Life;
									float wd = sqrt((float)(x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2))/dval3;
									if(wd < dminP)
									{
										dminP = wd;
									}
								}
								float dminN = 1.0e10; //a big number
								for(int i2=0; i2<vn.size(); ++i2)
								{
									int x3 = vn[i2].m_X;
									int y3 = vn[i2].m_Y;
									int z3 = vn[i2].m_Z;
									float dval3 = vn[i2].m_Life;
									float wd = sqrt((float)(x3-x2)*(x3-x2)+(y3-y2)*(y3-y2)+(z3-z2)*(z3-z2))/dval3;
									if(wd < dminN)
									{
										dminN = wd;
									}
								}
								if(dminP < dminN)
								{
									vc2.push_back(pn);
								}
							}
						}
					}
				}
			}
		}
		printf("Constrained-growth: #vc = %d, #vc2 = %d\n", vc.size(), vc2.size());
		vc = vc2;
	}
}

void
YetNewRegionPartition(vector<unsigned char>& C,
					  const vector<unsigned char>& B, //free-growth
					  const vector<float>& D,	//distance map computed on the free-growh
					  vector<unsigned char>& P,	//local maxima of D
					  int x, int y, int z,		//seed point
					  const int* dims)
{
	int i, j, k;
	vector<CParticle> vp;
	vector<CParticle> vn;

	float xf = x, yf = y, zf = z;
	float dthres = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);
	for(i=0; i<dims[2]; ++i)
	{
		for(j=0; j<dims[1]; ++j)
		{
			for(k=0; k<dims[0]; ++k)
			{
				if(GetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					float d = sqrt((float)(xf-k)*(xf-k)+(yf-j)*(yf-j)+(zf-i)*(zf-i));
					float dval = GetData3(D, k, j, i, dims[0], dims[1], dims[2], (float)0);
					if(d < dthres)
					{
						vp.push_back(CParticle(k, j, i, dval));
						SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)3);
					}
					else
					{
						vn.push_back(CParticle(k, j, i, dval));
						SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)2);
					}
				}
			}
		}
	}
	CParticle sd(x, y, z, dthres);
	if(find(vp.begin(), vp.end(), sd) == vp.end())
	{
		vp.insert(vp.begin(), sd);
	}

	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<dims[1]; ++j)
		{
			for(int k=0; k<dims[0]; ++k)
			{
				if(GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					float minval = -1;
					bool positive = false;
					for(int m=0; m<vp.size(); ++m)
					{
						int x1 = vp[m].m_X;
						int y1 = vp[m].m_Y;
						int z1 = vp[m].m_Z;
						float d = sqrt((float)(x1-k)*(x1-k)+(y1-j)*(y1-j)+(z1-i)*(z1-i));
						float cval = d/vp[m].m_Life;
						if(minval < 0 || cval < minval)
						{
							minval = cval;
							positive = true;
						}
					}
					for(int m=0; m<vn.size(); ++m)
					{
						int x1 = vn[m].m_X;
						int y1 = vn[m].m_Y;
						int z1 = vn[m].m_Z;
						float d = sqrt((float)(x1-k)*(x1-k)+(y1-j)*(y1-j)+(z1-i)*(z1-i));
						float cval = d/vn[m].m_Life;
						if(minval < 0 || cval < minval)
						{
							minval = cval;
							positive = false;
						}
					}
					if(positive)
					{
						SetData3(C, k, j, i, dims[0], dims[1], dims[2], (unsigned char)ForegroundColor);
					}
				}
			}
		}
	}
}

/*void
ConstrainedGrowth(vector<unsigned char>& B,
				  const vector<float>& D,
				  vector<unsigned char>& P,
				  int x, int y, int z,
				  const int* dims)
{
	const int PosLabel = 1; //label for a nodule
	const int NegLabel = 2; //label for a non-nodule
	float dthres = GetData3(D, x, y, z, dims[0], dims[1], dims[2], (float)0);

	int i, j, k;
	vector<CParticle> vpP(1, CParticle(x, y, z));
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

	incrementalTrace(B, D, vpP, dims);		
	incrementalTrace(Bn, D, vpN, dims);		

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
					if(mindsp > dthres && mindsp > mindsn)
					{
						SetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0);
					}
				}
			}
		}
	}
}*/

