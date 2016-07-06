#ifdef MEX_DLL
#include <mex.h>
#endif
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szConnectedComponent.h>

/*void
Initialize4SolidExp(const vector<unsigned short>& A,
					vector<float>& Q) 
{
  //double sgm2=2.*600.0*600.0;
  double sgm2=2.*500.0*500.0;
  double hv=1000.0;
  unsigned short usHv=1000;
  for(int i=0; i<A.size(); ++i) {
    if(A[i]>usHv)
      Q[i]=(float)1.0;
    else {
      double df=hv - (double)A[i];
      double s=exp(-df*df/sgm2);
      Q[i]=(float)s;
    }
  }
}*/

/*
extract everything that are reasonably bright.
*/
void
Initialize4NonBackground(const vector<unsigned short>& A,
					     vector<float>& Q) 
{
  unsigned short usHv=350;
  unsigned short usLv=150;
  float hv=(float) usHv;
  float lv=(float) usLv;
  for(int i=0; i<A.size(); ++i) {
    if(A[i]>usHv)
      Q[i]=(float)1.0;
    if(A[i]<usLv)
      Q[i]=(float)0.0;
    else {
      Q[i]=((float)A[i] - lv)/(hv - lv);
    }
  }
}

void
Initialize4SolidLinear(const vector<unsigned short>& A,
					   vector<float>& Q) 
{
  //double sgm2=2.*600.0*600.0;
  //double sgm2=2.*500.0*500.0;
  unsigned short usHv=1000;
  float hv=(float) usHv;
  for(int i=0; i<A.size(); ++i) {
    if(A[i]>usHv)
      Q[i]=(float)1.0;
    else {
      Q[i]=(float)A[i]/hv;
    }
  }
}

void
Initialize4SolidLinearLow(const vector<unsigned short>& A,
						  vector<float>& Q) 
{
  //double sgm2=2.*600.0*600.0;
  //double sgm2=2.*500.0*500.0;
  unsigned short usHv=800;
  float hv=(float) usHv;
  for(int i=0; i<A.size(); ++i) {
    if(A[i]>usHv)
      Q[i]=(float)1.0;
    else {
      Q[i]=(float)A[i]/hv;
    }
  }
}

void
Initialize4NonSolid(const vector<unsigned short>& A,
					vector<float>& Q) 
{
  float mv=400.0;
  float low=100; //200;
  float high=800; //1000; //900;
  float coeff = 0.2;
  for(int i=0; i<A.size(); ++i) {
    if(A[i]<mv)
      Q[i]=Max(0,((float)A[i]-low)/(mv-low));
    else 
      Q[i]=Max(0,(high-(float)A[i])/(high-mv));
  }
}

void
Initialize4NonSolidLow(const vector<unsigned short>& A,
					   vector<float>& Q) 
{
  float mv=300.0;
  float low=100; //200;
  float high=700; //1000; //900;
  float coeff = 0.2;
  for(int i=0; i<A.size(); ++i) {
    if(A[i]<mv)
      Q[i]=Max(0,((float)A[i]-low)/(mv-low));
    else 
      Q[i]=Max(0,(high-(float)A[i])/(high-mv));
	Q[i] = 0.5 + coeff * (Q[i] - 0.5);
  }
}

/*void
Initialize4PartSolid(const vector<unsigned short>& A,
					 vector<float>& Q) 
{
  float mv=300.0;
  float low=100; //200;
  float high=900; //900;
  float coeff = 0.2;
  for(int i=0; i<A.size(); ++i) {
    if(A[i]<mv)
      Q[i]=Max(0,((float)A[i]-low)/(mv-low));
    else 
      Q[i]=Max(0,(high-(float)A[i])/(high-mv));
	Q[i] = 0.5 + coeff * (Q[i] - 0.5);
  }
}*/

void
Initialize4PartSolid2(const vector<unsigned short>& A,
					  vector<float>& Q) {
  float mv=500.0;
  float low=100; //200;
  //float high=1000; //900;
  for(int i=0; i<A.size(); ++i) {
    if((float)A[i]<mv)
      Q[i]=Max(0,((float)A[i]-low)/(mv-low));
    else 
      Q[i]=1.0;
  }
}

void
ReactionDiffusion(vector<float>& A, 
                  int xD, 
                  int yD, 
                  int zD, 
                  float lambda,
				  float dfspeed[6],
				  float reacWeight) {
  int nump=xD*yD*zD;
  float v[6];
  float fzero=0;
  //float weight = 1.0;
  //printf("Using weight of %f\n", weight);
  int i, j, k;
  float dfspeed_sum = 0;
  for(i=0; i<6; ++i)
  {
	  dfspeed_sum += dfspeed[i];
  }
  for(k=0; k<zD; ++k) {
    for(j=0; j<yD; ++j) {
      for(i=0; i<xD; ++i) {
        float cval=GetData3(A,i,j,k,xD,yD,zD,fzero);
		float cval2 = 1.0 - cval;
        v[0]=GetData3(A,i-1,j,k,xD,yD,zD,cval);
        v[1]=GetData3(A,i+1,j,k,xD,yD,zD,cval);
        v[2]=GetData3(A,i,j-1,k,xD,yD,zD,cval);
        v[3]=GetData3(A,i,j+1,k,xD,yD,zD,cval);
        v[4]=GetData3(A,i,j,k-1,xD,yD,zD,cval);
        v[5]=GetData3(A,i,j,k+1,xD,yD,zD,cval);
        float difval1=0;
        float difval2=0;
        for(int m=0; m<6; ++m) {
          difval1+=dfspeed[m]*v[m];
          difval2+=dfspeed[m]*(1.-v[m]);
        }
        difval1/=dfspeed_sum;
        difval2/=dfspeed_sum;

		float sum = reacWeight*cval*cval + (1-cval)*(1-cval);
		float reacval1 = reacWeight*cval*cval/sum;
		float reacval2 = (1-cval)*(1-cval)/sum;
        float newval1=(1-lambda)*reacval1+lambda*difval1;
        float newval2=(1-lambda)*reacval2+lambda*difval2;
        float sum2=reacWeight*newval1*newval1+newval2*newval2; //this and the next may not be necessary
        float newval3=reacWeight*newval1*newval1/sum2;
        SetData3(A,i,j,k,xD,yD,zD,newval1);
      }
    }
  }
}

bool
FigureGroundSeparation(const vector<unsigned short>& A,  //input sub-volume
                       vector<unsigned char>& L,  //output segmentation
					   vector<float>& Q, //output probability image
                       int mode,                  //0: solid, 1: non-solid
                       int numiter,
                       float lambda,
                       const int* dimsA)
{   

  //mexPrintf("Casting short to float...\n");
  int i;

  //printf("Initializing the volume...\n");
  //initialize the state value
  float dfspeed[]={1,1,1,1,0.5,0.5};
  float reacWeight=1.0;
  if(mode==SZ_LABEL_SOLID_NODULE)
  {
    Initialize4SolidLinear(A, Q);
	lambda = 0.5; //default
  }
  else if(mode==SZ_LABEL_NONSOLID_NODULE)
  {
	Initialize4NonSolid(A, Q);
	lambda = 0.5; //default
	reacWeight = 1.0; //default
  }
  else if(mode==SZ_LABEL_PARTSOLID_NODULE)
  {
	Initialize4PartSolid2(A, Q);
	lambda = 0.5; //default
	reacWeight = 1.0; //default
  }
  else if(mode==SZ_LABEL_SOLID_NODULE_LOW)
  {
    Initialize4SolidLinearLow(A, Q);
	lambda = 0.5; //default
  }
  else if(mode==SZ_LABEL_NONSOLID_NODULE_LOW)
  {
	Initialize4NonSolidLow(A, Q);
	lambda = 0.5; //default
	reacWeight = 1.0; //default
  }
  else if(mode==SZ_LABEL_NON_BACKGROUND)
  {
	Initialize4NonBackground(A, Q);
	lambda = 0.5; //default
	reacWeight = 1.0; //default
  }
  else //assume solid
  {
    Initialize4SolidLinear(A, Q);
	lambda = 0.5;
	reacWeight = 1.0; //default
  }

  for(int iter=0; iter<numiter; ++iter) {
	  ReactionDiffusion(Q,dimsA[0],dimsA[1],dimsA[2],lambda, dfspeed, reacWeight);
  }

  for(i=0; i<Q.size(); ++i) {
	  if(Q[i]>.5) L[i]=ForegroundColor;
	  else L[i]=0;
  }

  return true; //right now, always returns true...
}

