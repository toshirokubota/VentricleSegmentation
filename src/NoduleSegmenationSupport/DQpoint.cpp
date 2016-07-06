//#include "stdafx.h"
#include <DQpoint.h>
#ifdef MEX_DLL
#include <mex.h>
#endif
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>

// Standard Output
ostream& operator<<(ostream& s, const DQpoint& v) 
{
	s << v.GetX() << " " << v.GetY() << " " << v.GetZ();
	s << " " << v.GetA() << " " << v.GetB() << " " << v.GetC();
	s << " " << v.GetD() << " " << v.GetE() << " " << v.GetF() <<"\n";
	return s;
}

istream& operator>>(istream& s, DQpoint& v) 
{
	real val;
	s >> val; v.SetX(val);
	s >> val; v.SetY(val);
	s >> val; v.SetZ(val);
	s >> val; v.SetA(val);
	s >> val; v.SetB(val);
	s >> val; v.SetC(val);
	s >> val; v.SetD(val);
	s >> val; v.SetE(val);
	s >> val; v.SetF(val);
	return s;
}

void
InitializeDensityQ(DQpointImage& res, const RealImage& image) 
{
	res = DQpointImage(image.NumBands(),image.NumRows(),image.NumCols());
	for(int k=0; k<image.NumBands(); ++k) 
	{
		for(int i=0; i<image.NumRows(); ++i) 
		{
			for(int j=0; j<image.NumCols(); ++j) 
			{
				real v00=image.GetPixel(k,i,j);
				DQpoint sp((real)j,(real)i,v00,0,0,0,0,0,0);
				res.SetPixel(k,i,j,sp);
			}
		}
	}
}

void
InitializeDensityQ(DQpointImage& res, const vector<real>& A, const int* dims) 
{
	res = DQpointImage(dims[2], dims[1], dims[0]);
	for(int k=0; k<res.NumBands(); ++k) 
	{
		for(int i=0; i<res.NumRows(); ++i) 
		{
			for(int j=0; j<res.NumCols(); ++j) 
			{
				real v00=GetData3(A, j, i, k, dims[0], dims[1], dims[2], (real)0);
				DQpoint sp((real)j,(real)i,v00,0,0,0,0,0,0);
				res.SetPixel(k,i,j,sp);
			}
		}
	}
}

void GetGammaMatrix(real G[MatDimDQ][MatDimDQ], real r, real t)
{
	G[0][0] = (t*t+4.0*r)/(t*t+8.0*r)/6.0;
	G[0][1] =-1/(t*t+8.0*r)/6.0;
	G[0][2] = 0.0;
	G[0][3] = -1/(t*t+8.0*r)/6.0;
	G[0][4] = 0.0;
	G[0][5] = -1/(t*t+8.0*r)/6.0;
	G[0][6] = 0.0;
	G[1][0] = -1/(t*t+8.0*r)/6.0;
	G[1][1] = 2.0/3.0*(t*t+6.0*r)/(t*t)/(t*t*t*t+12.0*r*t*t+32.0*r*r);
	G[1][2] = 0.0;
	G[1][3] = 1/(t*t*t*t+12.0*r*t*t+32.0*r*r)/6.0;
	G[1][4] = 0.0;
	G[1][5] = 1/(t*t*t*t+12.0*r*t*t+32.0*r*r)/6.0;
	G[1][6] = 0.0;
	G[2][0] = 0.0;
	G[2][1] = 0.0;
	G[2][2] = 1/(t*t+2.0*r)/2.0;
	G[2][3] = 0.0; G[2][4] = 0.0;
	G[2][5] = 0.0;
	G[2][6] = 0.0;
	G[3][0] = -1/(t*t+8.0*r)/6.0;
	G[3][1] = 1/(t*t*t*t+12.0*r*t*t+32.0*r*r)/6.0;
	G[3][2] = 0.0;
	G[3][3] = 2.0/3.0*(t*t+6.0*r)/(t*t)/(t*t*t*t+12.0*r*t*t+32.0*r*r);
	G[3][4] = 0.0;
	G[3][5] = 1/(t*t*t*t+12.0*r*t*t+32.0*r*r)/6.0;
	G[3][6] = 0.0;
	G[4][0] = 0.0;
	G[4][1] = 0.0;
	G[4][2] = 0.0;
	G[4][3] = 0.0;
	G[4][4] = 1/(t*t+2.0*r)/2.0;
	G[4][5] = 0.0;
	G[4][6] = 0.0;
	G[5][0] = -1/(t*t+8.0*r)/6.0;
	G[5][1] = 1/(t*t*t*t+12.0*r*t*t+32.0*r*r)/6.0;
	G[5][2] = 0.0;
	G[5][3] = 1/(t*t*t*t+12.0*r*t*t+32.0*r*r)/6.0;
	G[5][4] = 0.0;
	G[5][5] = 2.0/3.0*(t*t+6.0*r)/(t*t)/(t*t*t*t+12.0*r*t*t+32.0*r*r);
	G[5][6] = 0.0;
	G[6][0] = 0.0;
	G[6][1] = 0.0;
	G[6][2] = 0.0;
	G[6][3] = 0.0;
	G[6][4] = 0.0;
	G[6][5] = 0.0;
	G[6][6] = 1/(t*t+2.0*r)/2.0;
}

DQpoint
FindOptimumDensityQ(const DQpoint& fn, const DQpoint& fw, 
					const DQpoint& fe, const DQpoint& fs,
					const DQpoint& ff, const DQpoint& fb,
					real G[MatDimDQ][MatDimDQ], real H[3][3], 
					real x, real y) 
{
	real u[MatDimDQ] = {0, 0, 0, 0, 0, 0, 0};
	int i,j;

	//North-South - A, B
	{
		real x0[3]={fn.GetZ(),fn.GetA(),fn.GetB()};
		real x2[3]={fs.GetZ(),fs.GetA(),fs.GetB()};
		int id[3] = {0, 1, 2};
		for(i=0; i<3; ++i) {
			real sum=.0;
			for(j=0; j<3; ++j) {
				sum+=H[i][j]*x2[j]+H[j][i]*x0[j];
			}
			u[id[i]]+=sum;
		}

	}
	//West-East - C, D
	{
		real x0[3]={fw.GetZ(),fw.GetC(),fw.GetD()};
		real x2[3]={fe.GetZ(),fe.GetC(),fe.GetD()};
		int id[3] = {0, 3, 4};
		for(i=0; i<3; ++i) {
			real sum=.0;
			for(j=0; j<3; ++j) {
				sum+=H[i][j]*x2[j]+H[j][i]*x0[j];
			}
			u[id[i]]+=sum;
		}

	}
	//Back-Front - E, F
	{
		real x0[3]={fb.GetZ(),fb.GetE(),fb.GetF()};
		real x2[3]={ff.GetZ(),ff.GetE(),ff.GetF()};
		int id[3] = {0, 5, 6};
		for(i=0; i<3; ++i) {
			real sum=.0;
			for(j=0; j<3; ++j) {
				sum+=H[i][j]*x2[j]+H[j][i]*x0[j];
			}
			u[id[i]]+=sum;
		}

	}
	real b[MatDimDQ];
	for(i=0; i<MatDimDQ; ++i) {
		real sum=.0;
		for(j=0; j<MatDimDQ; ++j) {
			sum+=G[i][j]*u[j];
		}
		b[i]=-sum;
	}
	DQpoint fnew(x,y,b[0],b[1],b[2],b[3],b[4],b[5],b[6]);
	//DQpoint fnew(x,y,0,b[1],b[2],b[3],b[4],b[5],b[6]);
	return fnew;
}

void
UpdateQuadraticDensity(DQpointImage& density, 
					   real r, real t) 
{
	int i,j,k;
	real G[MatDimDQ][MatDimDQ]; //={{.5*a/b,-.5/b,0},{-.5/b,1./(t*t*b),0},{0,0,.5/c}};
	GetGammaMatrix(G, r, t);
	real H[3][3]={{-2.,-t*t,t},{-t*t,.0,-2.*r*t},{-t,2.*r*t,-2.*r}};
	for(k=1; k<density.NumBands()-1; ++k) 
	{
		for(i=1; i<density.NumRows()-1; ++i) 
		{
			for(j=1; j<density.NumCols()-1; ++j) 
			{
				DQpoint fi=density.GetPixel(k,i,j);
				DQpoint ff=density.GetPixelRepeat(k+1,i,j);
				DQpoint fb=density.GetPixelRepeat(k-1,i,j);
				DQpoint fn=density.GetPixelRepeat(k,i-1,j);
				DQpoint fw=density.GetPixelRepeat(k,i,j-1);
				DQpoint fe=density.GetPixelRepeat(k,i,j+1);
				DQpoint fs=density.GetPixelRepeat(k,i+1,j);
				DQpoint fnew=FindOptimumDensityQ(fn,fw, fe, fs, fb, ff,  
					G, H, (real)j, (real)i);
				density.SetPixel(k,i,j,fnew);
			}
		}
	}
}

void
UpdateQuadraticDensityMask(DQpointImage& density, 
						   real r, real t,
						   const vector<unsigned char>& M) 
{
	int i,j,k,n;
	real G[MatDimDQ][MatDimDQ]; //={{.5*a/b,-.5/b,0},{-.5/b,1./(t*t*b),0},{0,0,.5/c}};
	GetGammaMatrix(G, r, t);
	real H[3][3]={{-2.,-t*t,t},{-t*t,.0,-2.*r*t},{-t,2.*r*t,-2.*r}};
	const int dims[3] = {density.NumCols(), density.NumRows(), density.NumBands()};
	for(k=1; k<density.NumBands()-1; ++k) 
	{
		for(i=1; i<density.NumRows()-1; ++i) 
		{
			for(j=1; j<density.NumCols()-1; ++j) 
			{
				DQpoint fi=density.GetPixel(k,i,j);
				unsigned char mi = GetData3(M, j, i, k, dims[0], dims[1], dims[2], (unsigned char)0);
				DQpoint fn[NumNeighbors6];
				for(n=0; n<NumNeighbors6; ++n)
				{
					unsigned char mj = GetData3(M, j+XOffset6[n], i+YOffset6[n], k+ZOffset6[n], dims[0], dims[1], dims[2], (unsigned char)0);
					if(mi == mj)
					{
						fn[n] = density.GetPixelRepeat(k+ZOffset6[n], i+YOffset6[n], j+XOffset6[n]);
					}
					else
					{
						fn[n]  = DQpoint(j+XOffset6[n], i+YOffset[n], fi.GetZ(), 0, 0, 0, 0, 0, 0);
					}
				}
				DQpoint fnew=FindOptimumDensityQ(fn[0], fn[1], fn[2], fn[3], fn[4], fn[5],  
					G, H, (real)j, (real)i);
				density.SetPixel(k,i,j,fnew);
			}
		}
	}
}

