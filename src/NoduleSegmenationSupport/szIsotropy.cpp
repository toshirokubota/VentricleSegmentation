#ifdef CAD
#include "CADUtilityFunctions.h"
#endif

//#include <mex.h>

#include <math.h>
#include "szIsotropy.h"
#include "szMiscOperations.h"
#include "szMexUtilityTemplate.h"

void
IsotropyMeasure(vector<float>& D, 
                vector<unsigned char>& L,
                vector<float>& vWall, 
                vector<float>& vVessel, 
                vector<float>& vNodule, 
                const int* dims) {

  //printf("IsotropyMeasure\n");
  int i,j,k;
  int w=2;
  float alpha=5.0;
  bool fDescend=true;
  for(i=0; i<dims[2]; ++i) { 
    for(j=0; j<dims[1]; ++j) { 
      for(k=0; k<dims[0]; ++k) { 
        if(GetData3(L,k,j,i,dims[0],dims[1],dims[2],(unsigned char)0)) {
          vector<float> vst = computeStructuralTensor(D,L,k,j,i,w,dims);
          //vector<float> vst = computeHessian(D,L,k,j,i,w,dims);
          vector<float> ev(3);
          float egv1, egv2, egv3;
	        CAD::Eigen3D(vst[0], vst[1], vst[2], vst[4], vst[5], vst[8], egv1, egv2, egv3);
          ev[0]=Abs(egv1); ev[1]=Abs(egv2); ev[2]=Abs(egv3);
          sort(ev.begin(),ev.end());
          //printf("%d %d %d %f %f %f\n", k, j, i, egv1, egv2, egv3);
          //SetData3(vWall,k,j,i,dims[0],dims[1],dims[2],ev[0]);
          //SetData3(vVessel,k,j,i,dims[0],dims[1],dims[2],ev[1]);
          //SetData3(vNodule,k,j,i,dims[0],dims[1],dims[2],ev[2]);
          
          if(ev[2]>0) {
            /*float nev1=ev[1]/ev[2];
            float nev2=ev[0]/ev[2];
            //float pwall=exp(-alpha*(nev1+nev2));
            //float pvessel=exp(-alpha*((1-nev1)+nev2));
            //float pnodule=exp(-alpha*((1-nev1)+(1-nev2)));

            //transform the ratios non-linearly
            float t1 = alpha*(nev1*2-1.0);
            float t2 = alpha*(nev2*2-1.0);
            float pwall=1.0/(1+exp(+t1+t2));
            float pvessel=1.0/(1+exp(-t1+t2));
            float pnodule=1.0/(1+exp(-t1-t2));
            float sum=pwall+pvessel+pnodule;*/

            float beta=ev[2]/2;
            float beta2=5000;
            float dWall=ev[0]*ev[0]/beta2+ev[1]*ev[1]/beta2;
            float dNodule=(ev[0]-beta)*(ev[0]-beta)/(ev[2]*ev[2])+(ev[1]-beta)*(ev[1]-beta)/(ev[2]*ev[2]);
            float dVessel=ev[0]*ev[0]/beta2+(ev[1]-beta)*(ev[1]-beta)/(ev[2]*ev[2]);
            float sWall = exp(-dWall);
            float sNodule = exp(-dNodule);
            float sVessel = exp(-dVessel);
            float sum = sWall + sNodule + sVessel;
            
            SetData3(vWall,k,j,i,dims[0],dims[1],dims[2],sWall/sum);
            SetData3(vVessel,k,j,i,dims[0],dims[1],dims[2],sVessel/sum);
            SetData3(vNodule,k,j,i,dims[0],dims[1],dims[2],sNodule/sum);
          }
        }
      }
    }
  }
}

