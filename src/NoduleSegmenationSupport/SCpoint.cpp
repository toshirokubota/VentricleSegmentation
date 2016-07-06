//#include <stdafx.h>
#include <SCpoint.h>
#include <szMexUtility.h>

// Standard Output
ostream& operator<<(ostream& s, const SCpoint& v) {
  s << v.GetX() << " " << v.GetY() << " " << v.GetZ();
  s << " " << v.GetA() << " " << v.GetB() << " " << v.GetC();
  s << " " << v.GetD() << " " << v.GetE() << " " << v.GetF();
  s << " " << v.GetG() << " " << v.GetH() << " " << v.GetP() <<"\n";
  return s;
}

istream& operator>>(istream& s, SCpoint& v) {
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
  s >> val; v.SetG(val);
  s >> val; v.SetH(val);
  s >> val; v.SetP(val);
  return s;
}

SCpointImage
InitializeSurfaceC(const RealImage& image) {
   SCpointImage res(image.NumBands(),image.NumRows(),image.NumCols());
   for(int k=0; k<image.NumBands(); ++k) {
      for(int i=0; i<image.NumRows(); ++i) {
         for(int j=0; j<image.NumCols(); ++j) {
            real v00=image.GetPixel(k,i,j);
            SCpoint sp((real)j,(real)i,v00);
            res.SetPixel(k,i,j,sp);
         }
      }
   }
   return res;
}

/*
reads 10x10 Gamma matrices - which capture clique potential of the
quadratic spline energy model.
*/
int
ReadGammaMatricesC(real* g1,real* g2,real* g3,real* g4, 
                  real& tau, real& lambda,
                  char* filename) {
   ifstream in(filename,ios::in);
   if (!in){
      cout << "Error opening input file: " << filename << endl;
      return -1;
   }
   int i,j,dim=MatDimC*MatDimC;
   in >> tau;
   in >> lambda;
   cout << "Tau=" << tau << " Lambda=" << lambda << endl;
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

SCpoint
FindOptimumSurf(const SCpoint& fn, const SCpoint& fw, const SCpoint& fe, const SCpoint& fs,
                real* gn, real* gw, real* ge, real* gs, real x, real y) {
   real a[MatDimC];
   real b[MatDimC]={0,0,0,0,0,0,0,0,0,0};
   a[0]=fn.GetZ(); a[1]=fn.GetA(); a[2]=fn.GetB(); a[3]=fn.GetC();
   a[4]=fn.GetD(); a[5]=fn.GetE(); a[6]=fn.GetF(); a[7]=fn.GetG();
   a[8]=fn.GetH(); a[9]=fn.GetP();
   int i,j;
   real w=-1.;
   for(i=0; i<MatDimC; ++i) {
      for(j=0; j<MatDimC; ++j) {
         b[i]+=w*gn[MatDimC*i+j]*a[j];
      }
   }
   a[0]=fw.GetZ(); a[1]=fw.GetA(); a[2]=fw.GetB(); a[3]=fw.GetC();
   a[4]=fw.GetD(); a[5]=fw.GetE(); a[6]=fw.GetF(); a[7]=fw.GetG();
   a[8]=fw.GetH(); a[9]=fw.GetP();
   for(i=0; i<MatDimC; ++i) {
      for(j=0; j<MatDimC; ++j) {
         b[i]+=w*gw[MatDimC*i+j]*a[j];
      }
   }
   a[0]=fe.GetZ(); a[1]=fe.GetA(); a[2]=fe.GetB(); a[3]=fe.GetC();
   a[4]=fe.GetD(); a[5]=fe.GetE(); a[6]=fe.GetF(); a[7]=fe.GetG();
   a[8]=fe.GetH(); a[9]=fe.GetP();
   for(i=0; i<MatDimC; ++i) {
      for(j=0; j<MatDimC; ++j) {
         b[i]+=w*ge[MatDimC*i+j]*a[j];
      }
   }
   a[0]=fs.GetZ(); a[1]=fs.GetA(); a[2]=fs.GetB(); a[3]=fs.GetC();
   a[4]=fs.GetD(); a[5]=fs.GetE(); a[6]=fs.GetF(); a[7]=fs.GetG();
   a[8]=fs.GetH(); a[9]=fs.GetP();
   for(i=0; i<MatDimC; ++i) {
      for(j=0; j<MatDimC; ++j) {
         b[i]+=w*gs[MatDimC*i+j]*a[j];
      }
   }
   SCpoint fnew(x,y,b[0],b[1],b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9]);
   return fnew;
}


/*
Check the discrepancy of f2 from f1.
f=Z+A*s^2+B*s+C*t^2+D*t+E*s*t
*/
real
CompatibilityMeasure(const SCpoint& f1, const SCpoint& f2,
                     real s1, real t1, real s2, real t2, real lambda) {
   real z1=f1.Evaluate(s1,t1);
   real z2=f2.Evaluate(s2,t2);
   Complex td=f1.Tangent(s1,t1)-f2.Tangent(s2,t2);
   real d2=td.Power();
   return (z1-z2)*(z1-z2)+lambda*d2;
}

void
UpdateCubicSurface(SCpointImage& surf,
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda) {
   int i,j,k;
   for(k=0; k<surf.NumBands(); ++k) {
      for(i=1; i<surf.NumRows()-1; ++i) {
         for(j=1; j<surf.NumCols()-1; ++j) {
            SCpoint fi=surf.GetPixel(k,i,j);
            SCpoint fn=surf.GetPixelRepeat(k,i-1,j);
            SCpoint fw=surf.GetPixelRepeat(k,i,j-1);
            SCpoint fe=surf.GetPixelRepeat(k,i,j+1);
            SCpoint fs=surf.GetPixelRepeat(k,i+1,j);
            SCpoint fnew=FindOptimumSurf(fn,fw,fe,fs,gn,gw,ge,gs,(real)j,(real)i);
            surf.SetPixel(k,i,j,fnew);
         }
      }
   }
}

void
UpdateCubicSurfaceBreak(SCpointImage& surf,
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda, real thres) {
   int i,j,k;
   for(k=0; k<surf.NumBands(); ++k) {
      for(i=1; i<surf.NumRows()-1; ++i) {
         for(j=1; j<surf.NumCols()-1; ++j) {
            SCpoint fi=surf.GetPixel(k,i,j);
            SCpoint fn=surf.GetPixelRepeat(k,i-1,j);
			if(Abs(fi.GetZ() - fn.GetZ()) > thres)
			{
				fn = fi;
			}
            SCpoint fw=surf.GetPixelRepeat(k,i,j-1);
			if(Abs(fi.GetZ() - fw.GetZ()) > thres)
			{
				fw = fi;
			}
            SCpoint fe=surf.GetPixelRepeat(k,i,j+1);
			if(Abs(fi.GetZ() - fe.GetZ()) > thres)
			{
				fe = fi;
			}
            SCpoint fs=surf.GetPixelRepeat(k,i+1,j);
			if(Abs(fi.GetZ() - fs.GetZ()) > thres)
			{
				fs = fi;
			}
            SCpoint fnew=FindOptimumSurf(fn,fw,fe,fs,gn,gw,ge,gs,(real)j,(real)i);
            surf.SetPixel(k,i,j,fnew);
         }
      }
   }
}

void
UpdateCubicSurfaceMask(SCpointImage& surf,
                        real* gn, real* gw, real* ge, real* gs,
                        real tau, real lambda, 
						const ByteImage& mask)
{
   int i,j,k;
   real total_sum = 0;
   for(k=0; k<surf.NumBands(); ++k) {
      for(i=1; i<surf.NumRows()-1; ++i) {
         for(j=1; j<surf.NumCols()-1; ++j) {
            SCpoint fi=surf.GetPixel(k,i,j);
			unsigned char mi = mask.GetPixel(k, i, j);
            SCpoint fn=surf.GetPixelRepeat(k,i-1,j);
			unsigned char mn = mask.GetPixelRepeat(k, i-1, j);
			if(mi!=mn)
			{
				fn = fi;
			}
            SCpoint fw=surf.GetPixelRepeat(k,i,j-1);
			unsigned char mw = mask.GetPixelRepeat(k, i, j-1);
			if(mi!=mw)
			{
				fw = fi;
			}
            SCpoint fe=surf.GetPixelRepeat(k,i,j+1);
			unsigned char me = mask.GetPixelRepeat(k, i, j+1);
			if(mi != me)
			{
				fe = fi;
			}
            SCpoint fs=surf.GetPixelRepeat(k,i+1,j);
			unsigned char ms = mask.GetPixelRepeat(k, i+1, j);
			if(mi != ms)
			{
				fs = fi;
			}
            SCpoint fnew=FindOptimumSurf(fn,fw,fe,fs,gn,gw,ge,gs,(real)j,(real)i);
            surf.SetPixel(k,i,j,fnew);

			/*real sum=0;
            sum+=CompatibilityMeasure(fnew,fn,0,-tau,0,0,lambda);
            sum+=CompatibilityMeasure(fnew,fn,0,0.,0,tau,lambda);
            sum+=CompatibilityMeasure(fnew,fw,-tau,0,0,0,lambda);
            sum+=CompatibilityMeasure(fnew,fw,0,0.,tau,0,lambda);
            sum+=CompatibilityMeasure(fnew,fe,tau,0,0,0,lambda);
            sum+=CompatibilityMeasure(fnew,fe,0,0.,-tau,0,lambda);
            sum+=CompatibilityMeasure(fnew,fs,0,tau,0,0,lambda);
            sum+=CompatibilityMeasure(fnew,fs,0,0.,0,-tau,lambda);
            total_sum+=sum;*/
         }
      }
   }
   //printf("total support: %f\n", total_sum);
}


CCurvature 
SCpoint::Curvature(real s, real t)
{
	Complex ct = Tangent(s, t);

	real gx=ct.GetReal();
	real gy=ct.GetImag();
	real gxx=6*GetA()*s + 2*GetB()*t + 2*GetE();;
	real gyy=6*GetD()*t + 2*GetC()*s + 2*GetG();;
	real gxy=2*GetB()*s + 2*GetC()*t + GetF();

	real K=(gxx*gyy-gxy*gxy)/pow(1.+gx*gx+gy*gy, 2.0);
	real H=.5*((1.+gx*gx)*gyy-2.*gx*gy*gxy+(1.+gy*gy)*gxx)/pow(1.+gx*gx+gy*gy,1.5);

	CCurvature cv;
	cv.K = K;
	cv.H = H;
	cv.k1 = H - sqrt(H*H-K);
	cv.k2 = H + sqrt(H*H-K);

	return cv;
}
