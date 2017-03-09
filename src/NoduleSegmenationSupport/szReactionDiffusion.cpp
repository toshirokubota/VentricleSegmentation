#ifdef MEX_DLL
#include <mex.h>
#endif
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szConnectedComponent.h>


void
ReactionDiffusion4D(vector<float>& A,
	int xD,
	int yD,
	int zD,
	int tD,
	vector<float>& C,
	float dfspeed[8],
	float rate)
{
	int nump = xD*yD*zD*tD;
	float v[8];
	for (int m = 0; m < zD; ++m)
	{
		for (int k = 0; k < zD; ++k)
		{
			for (int j = 0; j < yD; ++j)
			{
				for (int i = 0; i < xD; ++i)
				{
					float cval = GetData4(A, i, j, k, m, xD, yD, zD, tD, 0.0f);
					v[0] = GetData4(A, i - 1, j, k, m, xD, yD, zD, tD, cval);
					v[1] = GetData4(A, i + 1, j, k, m, xD, yD, zD, tD, cval);
					v[2] = GetData4(A, i, j - 1, k, m, xD, yD, zD, tD, cval);
					v[3] = GetData4(A, i, j + 1, k, m, xD, yD, zD, tD, cval);
					v[4] = GetData4(A, i, j, k - 1, m, xD, yD, zD, tD, cval);
					v[5] = GetData4(A, i, j, k + 1, m, xD, yD, zD, tD, cval);
					v[6] = GetData4(A, i, j, k, m - 1, xD, yD, zD, tD, cval);
					v[7] = GetData4(A, i, j, k, m + 1, xD, yD, zD, tD, cval);
					float difval = 0;
					for (int m = 0; m < 8; ++m) {
						difval += dfspeed[m] * (v[m] - cval);
					}
					//difval/=dfspeed_sum;

					float fitval1 = C[0] * cval + C[1] * (1.0 - cval);
					float fitval2 = C[2] * cval + C[3] * (1.0 - cval);
					float avval = cval * fitval1 + (1.0 - cval) * fitval2;
					float reacval = cval * (fitval1 - avval);
					float dval = reacval + difval;
					float newval = cval + rate * dval;
					newval = Max(0.0, Min(1.0, newval));
					SetData4(A, i, j, k, m, xD, yD, zD, tD, newval);
					/*float dif = newval - cval;
					if(Abs(dif) > maxDif)
					{
					maxDif = Abs(dif);
					}*/
				}
			}
		}
	}
	//printf("ReactionDiffusionNew: max difference = %f\n", maxDif);
}



void
ReactionDiffusion3D(vector<float>& A, 
					 int xD, 
					 int yD, 
					 int zD,
					 vector<float>& C,
					 float dfspeed[6],
					 float rate) 
{
	int nump=xD*yD*zD;
	float v[6];
	float fzero=0;
	int i, j, k;
	float dfspeed_sum = 0;
	for(i=0; i<6; ++i)
	{
		dfspeed_sum += dfspeed[i];
	}
	//float maxDif = 0;
	for(k=0; k<zD; ++k) 
	{
		for(j=0; j<yD; ++j) 
		{
			for(i=0; i<xD; ++i) 
			{
				float cval=GetData3(A,i,j,k,xD,yD,zD,fzero);        
				v[0]=GetData3(A,i-1,j,k,xD,yD,zD,cval);
				v[1]=GetData3(A,i+1,j,k,xD,yD,zD,cval);
				v[2]=GetData3(A,i,j-1,k,xD,yD,zD,cval);
				v[3]=GetData3(A,i,j+1,k,xD,yD,zD,cval);
				v[4]=GetData3(A,i,j,k-1,xD,yD,zD,cval);
				v[5]=GetData3(A,i,j,k+1,xD,yD,zD,cval);
				float difval=0;
				for(int m=0; m<6; ++m) {
					difval+=dfspeed[m]*(v[m]-cval);
				}
				//difval/=dfspeed_sum;

				float fitval1 = C[0] * cval + C[1] * (1.0-cval);
				float fitval2 = C[2] * cval + C[3] * (1.0-cval);
				float avval = cval * fitval1 + (1.0-cval) * fitval2;
				float reacval = cval * (fitval1 - avval);
				float dval = reacval + difval;
				float newval = cval + rate * dval;
				newval = Max(0.0, Min(1.0, newval));
				SetData3(A,i,j,k,xD,yD,zD,newval);
				/*float dif = newval - cval;
				if(Abs(dif) > maxDif)
				{
					maxDif = Abs(dif);
				}*/
			}
		}
	}
	//printf("ReactionDiffusionNew: max difference = %f\n", maxDif);
}


