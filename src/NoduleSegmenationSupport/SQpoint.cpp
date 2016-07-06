//#include "stdafx.h"
#include <SQpoint.h> 
#include <szMexUtility.h>

CCurvature
SQpoint:: Curvature(real s, real t) const
{
	Complex ct = Tangent(s, t);

	real gx=ct.GetReal();
	real gy=ct.GetImag();
	real gxx=2*GetA();
	real gyy=2*GetC();
	real gxy=2*GetE();

	real K=(gxx*gyy-gxy*gxy)/pow(1.+gx*gx+gy*gy, 2.0);
	real H=.5*((1.+gx*gx)*gyy-2.*gx*gy*gxy+(1.+gy*gy)*gxx)/pow(1.+gx*gx+gy*gy,1.5);

	CCurvature cv;
	cv.K = K;
	cv.H = H;
	cv.k1 = H - sqrt(H*H-K);
	cv.k2 = H + sqrt(H*H-K);

	return cv;
}

// Standard Output
ostream& operator<<(ostream& s, const SQpoint& v) 
{
  s << v.GetX() << " " << v.GetY() << " " << v.GetZ();
  s << " " << v.GetA() << " " << v.GetB();
  s << " " << v.GetC() << " " << v.GetD() << " " << v.GetE() <<"\n";
  return s;
}

istream& operator>>(istream& s, SQpoint& v) 
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
  return s;
}

SQpointImage
InitializeSurfaceQ(const RealImage& image) 
{
   SQpointImage res(1,image.NumRows(),image.NumCols());
   for(int i=0; i<image.NumRows(); ++i) 
   {
      for(int j=0; j<image.NumCols(); ++j) 
	  {
         real v00=image.GetPixel(0,i,j);
         SQpoint sp((real)j,(real)i,v00,0,0,0,0,0);
         res.SetPixel(0,i,j,sp);
      }
   }
   return res;
}

int
ReadGammaMatricesQ(real* g1,real* g2,real* g3,real* g4, 
                  real& tau, real& lambda,
                  char* filename) 
{
   ifstream in(filename,ios::in);
   if (!in)
   {
      cout << "Error opening input file: " << filename << endl;
      return -1;
   }
   int i,j,dim=MatDimQ*MatDimQ;

   in>>tau;
   in>>lambda;
   for(i=0; i<dim; ++i) 
      in>>g1[i];
   
   for(i=0; i<dim; ++i) 
      in>>g2[i];

   for(i=0; i<dim; ++i) 
      in>>g3[i];

   for(i=0; i<dim; ++i) 
      in>>g4[i];

   return 0;
}

SQpoint
FindOptimumSurfQ(const SQpoint& fn, const SQpoint& fw, const SQpoint& fe, const SQpoint& fs,
                real* gn, real* gw, real* ge, real* gs, real x, real y) 
{
   real a[MatDimQ];
   real b[MatDimQ]={0,0,0,0,0,0};
   a[0]=fn.GetZ(); a[1]=fn.GetA(); a[2]=fn.GetB(); a[3]=fn.GetC();
   a[4]=fn.GetD(); a[5]=fn.GetE();
   int i,j;
   real w=-1.;
   for(i=0; i<MatDimQ; ++i) 
   {
      for(j=0; j<MatDimQ; ++j) 
	  {
         b[i]+=w*gn[6*i+j]*a[j];
      }
   }
   a[0]=fw.GetZ(); a[1]=fw.GetA(); a[2]=fw.GetB(); a[3]=fw.GetC();
   a[4]=fw.GetD(); a[5]=fw.GetE();
   for(i=0; i<MatDimQ; ++i) 
   {
      for(j=0; j<MatDimQ; ++j) 
	  {
         b[i]+=w*gw[MatDimQ*i+j]*a[j];
      }
   }
   a[0]=fe.GetZ(); a[1]=fe.GetA(); a[2]=fe.GetB(); a[3]=fe.GetC();
   a[4]=fe.GetD(); a[5]=fe.GetE();
   for(i=0; i<MatDimQ; ++i) 
   {
      for(j=0; j<MatDimQ; ++j) 
	  {
         b[i]+=w*ge[MatDimQ*i+j]*a[j];
      }
   }
   a[0]=fs.GetZ(); a[1]=fs.GetA(); a[2]=fs.GetB(); a[3]=fs.GetC();
   a[4]=fs.GetD(); a[5]=fs.GetE();
   for(i=0; i<MatDimQ; ++i) 
   {
      for(j=0; j<MatDimQ; ++j) 
	  {
         b[i]+=w*gs[MatDimQ*i+j]*a[j];
      }
   }
   SQpoint fnew(x,y,b[0],b[1],b[2],b[3],b[4],b[5]);
   return fnew;
}

void
UpdateQuadraticSurface(SQpointImage& surf, 
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda) 
{
   int i,j,k;
   for(i=1; i<surf.NumRows()-1; ++i) 
   {
      for(j=1; j<surf.NumCols()-1; ++j) 
	  {
         SQpoint fi=surf.GetPixel(0,i,j);
         SQpoint fn=surf.GetPixelRepeat(0,i-1,j);
         SQpoint fw=surf.GetPixelRepeat(0,i,j-1);
         SQpoint fe=surf.GetPixelRepeat(0,i,j+1);
         SQpoint fs=surf.GetPixelRepeat(0,i+1,j);
         SQpoint fnew=FindOptimumSurfQ(fn,fw,fe,fs,gn,gw,ge,gs,(real)j,(real)i);
         surf.SetPixel(0,i,j,fnew);
      }
   }
}

void
UpdateQuadraticSurface(SQpointImage& surf, 
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda,
						const ByteImage& mask) 
{
	int xoff[]={0, -1, 1, 0};
	int yoff[]={-1, 0, 0, 1};
   int i,j,k;
   for(i=0; i<surf.NumRows(); ++i) 
   {
      for(j=0; j<surf.NumCols(); ++j) 
	  {
         SQpoint fi=surf.GetPixel(0,i,j);
		 unsigned char mval = mask.GetPixel(0, i, j);
		 SQpoint fn[4];
		 for(int n=0; n<4; ++n)
		 {
			 unsigned char mval2 = mask.GetPixelDefault(0, i+yoff[n], j+xoff[n], mval);
			 if(mval == mval2)
			 {
				 fn[n] = surf.GetPixelRepeat(0, i+yoff[n], j+xoff[n]);
			 }
			 else
			 {
				 fn[n] = fi;
			 }
		 }
         SQpoint fnew=FindOptimumSurfQ(fn[0], fn[1], fn[2], fn[3],gn,gw,ge,gs,(real)j,(real)i);
         surf.SetPixel(0,i,j,fnew);
      }
   }
}

void
UpdateQuadraticSurfaceBreak(SQpointImage& surf,
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda, real thres) {
   int i,j,k;
   for(k=0; k<surf.NumBands(); ++k) {
      for(i=1; i<surf.NumRows()-1; ++i) {
         for(j=1; j<surf.NumCols()-1; ++j) {
            SQpoint fi=surf.GetPixel(k,i,j);
            SQpoint fn=surf.GetPixelRepeat(k,i-1,j);
			if(Abs(fi.GetZ() - fn.GetZ()) > thres)
			{
				fn = fi;
			}
            SQpoint fw=surf.GetPixelRepeat(k,i,j-1);
			if(Abs(fi.GetZ() - fw.GetZ()) > thres)
			{
				fw = fi;
			}
            SQpoint fe=surf.GetPixelRepeat(k,i,j+1);
			if(Abs(fi.GetZ() - fe.GetZ()) > thres)
			{
				fe = fi;
			}
            SQpoint fs=surf.GetPixelRepeat(k,i+1,j);
			if(Abs(fi.GetZ() - fs.GetZ()) > thres)
			{
				fs = fi;
			}
            SQpoint fnew=FindOptimumSurfQ(fn,fw,fe,fs,gn,gw,ge,gs,(real)j,(real)i);
            surf.SetPixel(k,i,j,fnew);
         }
      }
   }
}


