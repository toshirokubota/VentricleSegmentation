#include <smoothen.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
#include <SCpoint.h>
#include <DQpoint.h>

void
smoothen3(vector<float>& A, 
		  const vector<unsigned char>& M,
		  int niter,
		  real* gn, real* gs, real* gw, real* ge,
		  real tau, real lambda,
		  const int* dims)
{
	real maxDiff = 20;
	int i, j, k;
	//X-Y Planes
	{
		//printf("X-Y Planes:\n");
		RealImage image(1, dims[1], dims[0]);
		ByteImage mask(1, dims[1], dims[0]);
		for(k=0; k<dims[2]; ++k)
		{
			//printf("%d slice:\n", k);
			for(i=0; i<dims[1]; ++i)
			{
				for(j=0; j<dims[0]; ++j)
				{
					float val = GetData3(A, j, i, k, dims[0], dims[1], dims[2], (float)0);
					unsigned char mval = GetData3(M, j, i, k, dims[0], dims[1], dims[2], (unsigned char)0);
					mask.SetPixel(0, i, j, mval);
					val = Max(100, Min(1200, val));
					image.SetPixel(0, i, j, (double(val)));
				}
			}
			SCpointImage surf = InitializeSurfaceC(image);
			for(int m=0; m<niter; ++m)
			{
				UpdateCubicSurfaceMask(surf, gn, gw, ge, gs, tau, lambda, mask);
			}
			for(i=0; i<dims[1]; ++i)
			{
				for(j=0; j<dims[0]; ++j)
				{
					SetData3(A, j, i, k, dims[0], dims[1], dims[2], (float)surf.GetPixel(0, i, j).GetZ());
				}
			}
		}
	}
	//X-Z Planes
	{
		//printf("X-Z Planes:\n");
		RealImage image(1, dims[2], dims[0]);
		ByteImage mask(1, dims[1], dims[0]);
		for(k=0; k<dims[1]; ++k)
		{
			//printf("%d slice:\n", k);
			for(i=0; i<dims[2]; ++i)
			{
				for(j=0; j<dims[0]; ++j)
				{
					float val = GetData3(A, j, k, i, dims[0], dims[1], dims[2], (float)0);
					unsigned char mval = GetData3(M, j, k, i, dims[0], dims[1], dims[2], (unsigned char)0);
					mask.SetPixel(0, i, j, mval);
					val = Max(1000, Min(1200, val));
					image.SetPixel(0, i, j, (double(val)));
				}
			}
			SCpointImage surf = InitializeSurfaceC(image);
			for(int m=0; m<niter; ++m)
			{
				UpdateCubicSurfaceMask(surf, gn, gw, ge, gs, tau, lambda, mask);
			}
			for(i=0; i<dims[2]; ++i)
			{
				for(j=0; j<dims[0]; ++j)
				{
					SetData3(A, j, k, i, dims[0], dims[1], dims[2], (float)surf.GetPixel(0, i, j).GetZ());
				}
			}
		}
	}
	//Y-Z Planes
	/*{
		RealImage image(1, dims[2], dims[1]);
		ByteImage mask(1, dims[1], dims[0]);
		for(k=0; k<dims[0]; ++k)
		{
			for(i=0; i<dims[2]; ++i)
			{
				for(j=0; j<dims[1]; ++j)
				{
					float val = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (float)0);
					unsigned char mval = GetData3(M, j, i, k, dims[0], dims[1], dims[2], (unsigned char)0);
					mask.SetPixel(0, i, j, mval);
					val = Max(1000, Min(1200, val));
					image.SetPixel(0, i, j, (double(val)));
				}
			}
			SCpointImage surf = InitializeSurfaceC(image);
			for(int m=0; m<niter; ++m)
			{
				UpdateCubicSurfaceMask(surf, gn, gw, ge, gs, tau, lambda, mask);
			}
			for(i=0; i<dims[2]; ++i)
			{
				for(j=0; j<dims[1]; ++j)
				{
					SetData3(A, k, j, i, dims[0], dims[1], dims[2], (float)surf.GetPixel(0, i, j).GetZ());
				}
			}
		}
	}*/
}


void
smoothen3D(vector<float>& A, 
		  const vector<unsigned char>& M,
		  int niter,
		  real tau, real lambda,
		  const int* dims)
{

	DQpointImage density;

	InitializeDensityQ(density, A, dims);

	for(int n=0; n<niter; ++n)
	{
		//UpdateQuadraticDensityMask(density, lambda, tau, M); 
		UpdateQuadraticDensity(density, lambda, tau); 
	}
	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<dims[1]; ++j)
		{
			for(int k=0; k<dims[0]; ++k)
			{
				SetData3(A, k, j, i, dims[0], dims[1], dims[2], density.GetPixel(i, j, k).GetZ());
			}
		}
	}
}

void
smoothen3L(vector<float>& A, 
		   const vector<unsigned char>& M,
		   int niter,
		   float alpha,
		   const int* dims)
{
	for(int n=0; n<niter; ++n)
	{
		for(int i=0; i<dims[2]; ++i)
		{
			for(int j=0; j<dims[1]; ++j)
			{
				for(int k=0; k<dims[0]; ++k)
				{
					float sum = 0;
					float val0 = GetData3(A, k, j, i, dims[0], dims[1], dims[2], (float)0);
					for(int m=0; m<NumNeighbors6; ++m)
					{
						float val = GetData3(A, k+XOffset6[m], j+YOffset6[m], i+ZOffset6[m], dims[0], dims[1], dims[2], val0);
						sum += val;
					}
					float newval = alpha*val0 + (1-alpha)*sum/NumNeighbors6;

					SetData3(A, k, j, i, dims[0], dims[1], dims[2], (float)newval);
				}
			}
		}
	}
}

