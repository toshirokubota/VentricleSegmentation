//#include <mex.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szMyNeighborOp.h>
//#include <DistanceTransform.h>
//#include <LocalExtrema.h>
#include <szConnectedComponent.h>

/*
The case is degenerate if it is zero everywhere or non-zero everywhere
The routine assumes that L contains either 0 or 1.
*/
bool
IsDegenerate(const vector<unsigned char>& L) {
  for(int i=1; i<L.size(); ++i) {
    if(L[i]!=L[0])
      return false;
  }
  return true;
}

/*
The foreground is considered solitary, if 
1) the foreground volume is of significant size.
2) it does not touch any of six volume faces.
3) the fogreground is compact on its own
If so, we consider it to be non-solitary.
*/
bool
IsSolitary(const vector<unsigned char>& L, 
           int max_size, 
           const int* dims) {
  int i,j,k;
  //if the foreground is big, then it is highly likely that it is not solitary
  int sum=0;
  for(i=0; i<L.size(); ++i) {
    sum+=(L[i]>0? 1: 0);
  }
  if(sum>max_size) {
    return false;
  }
  //check if the foreground touches the VOI boundary
  int xD=dims[0];
  int yD=dims[1];
  int zD=dims[2];
  for(j=0; j<yD; ++j) {
    for(i=0; i<xD; ++i) {
      unsigned char cval=GetData3(L,i,j,0,dims[0],dims[1],dims[2],(unsigned char)0);
      if(cval>0) {
        return false;
      }
    }
  }
  for(j=0; j<yD; ++j) {
    for(i=0; i<xD; ++i) {
      unsigned char cval=GetData3(L,i,j,zD-1,dims[0],dims[1],dims[2],(unsigned char)0);
      if(cval>0) {
        return false;
      }
    }
  }
  for(j=0; j<zD; ++j) {
    for(i=0; i<xD; ++i) {
      unsigned char cval=GetData3(L,i,0,j,dims[0],dims[1],dims[2],(unsigned char)0);
      if(cval>0) {
        return false;
      }
    }
  }
  for(j=0; j<zD; ++j) {
    for(i=0; i<xD; ++i) {
      unsigned char cval=GetData3(L,i,yD-1,j,dims[0],dims[1],dims[2],(unsigned char)0);
      if(cval>0) {
        return false;
      }
    }
  }
  for(j=0; j<zD; ++j) {
    for(i=0; i<yD; ++i) {
      unsigned char cval=GetData3(L,0,i,j,dims[0],dims[1],dims[2],(unsigned char)0);
      if(cval>0) {
        return false;
      }
    }
  }
  for(j=0; j<zD; ++j) {
    for(i=0; i<yD; ++i) {
      unsigned char cval=GetData3(L,xD-1,i,j,dims[0],dims[1],dims[2],(unsigned char)0);
      if(cval>0) {
        return false;
      }
    }
  }

  //check if it is compact
  int minx=xD, miny=yD, minz=zD;
  int maxx=0, maxy=0, maxz=0;
  int count=0; //the number of foreground voxels
  for(j=0; j<zD; ++j) {
    for(i=0; i<yD; ++i) {
      for(k=0; k<xD; ++k) {
        unsigned char cval=GetData3(L,k,i,j,dims[0],dims[1],dims[2],(unsigned char)0);
        if(cval) {
          if(k<minx) minx=k;
          if(k>maxx) maxx=k;
          if(i<miny) miny=i;
          if(i>maxy) maxy=i;
          if(j<minz) minz=j;
          if(j>maxz) maxz=j;
          count++;
        }
      }
    }
  }

  //printf("%d %d %d %d %d %d %d\n", count, minx, maxx, miny, maxy, minz, maxz);
  int bbFG=(maxx-minx+1)*(maxy-miny+1)*(maxz-minz+1); //bounding box of the foreground
  if ((float)count/bbFG<.25)
    return false;

  return true;
}

/*
This is a wrapper for IsSolitary.  We check if a component that contains
a point at (x, y, z) is solitary or not.
*/
bool
IsSolitary(const vector<unsigned char>& L, 
           int x, int y, int z,
           int max_size, 
           const int* dims) 
{
  vector<int> C(L.size(), 0);
  int numClusters = ConnectedComponentAnalysisBigger(C, L, NeighborhoodFour, (unsigned char)0, 3, dims);
  int label = GetData3(C, x, y, z, dims[0], dims[1], dims[2], 0);
  //printf("num clusters = %d, label = %d\n", numClusters, label);
  if(label)
  {
    vector<unsigned char> S(L.size(), 0);
    for(int j=0; j<C.size(); ++j)
    {
      if(C[j] == label)
      {
        S[j] = 1;
      }
      else
      {
        S[j] = 0;
      }
    }
    if(IsSolitary(S, max_size, dims))
    {
      return true;
    }
    else
    {
      return false;
    }
  }
  else
  {
    return false;
  }
}

bool
IsWallAttached(const vector<unsigned char>& L,
               const int* dims) {
  vector<unsigned char> M=L;
  int i,j,k;
  //clear the inside the volume but leave the bounday intact
  for(k=1; k<dims[2]-1; ++k) {
    for(j=1; j<dims[1]-1; ++j) {
      for(i=1; i<dims[0]-1; ++i) {
        SetData3(M,i,j,k,dims[0],dims[1],dims[2],(unsigned char)0);
      }
    }
  }

  unsigned char value = 0;
  int ndim = 3;
  vector<int> C(L.size(),0);
  int numclusters = ConnectedComponentAnalysisBigger(C,M,NeighborhoodFour,value, ndim, dims);

  vector<int> vcount(numclusters,0);

  for(i=0; i<C.size(); ++i) {
    if(C[i])
      vcount[C[i]-1]++;
  }
  //now look through the counter 
  //It is wall attached if the biggest cluster is larger than the half the axial slice size
  for(i=0; i<numclusters; ++i) {
    if(vcount[i]>dims[0]*dims[1]/2)
      return true;
  }

  return false;

  //now check 8 corneres of the volume.
  //if wall-attached, a big cluster streatches over one of corners
  /*int xcorner[8]={0,dims[0]-1,0,dims[0]-1,0,dims[0]-1,0,dims[0]-1};
  int ycorner[8]={0,0,dims[1]-1,dims[1]-1,0,0,dims[1]-1,dims[1]-1};
  int zcorner[8]={0,0,0,0,dims[2]-1,dims[2]-1,dims[2]-1,dims[2]-1};
  vector<int> labels;

  for(int m=0; m<8; ++m) {
    int lb = GetData3(C,xcorner[m],ycorner[m],zcorner[m],dims[0],dims[1],dims[2],(int)0);
    if(lb) {
      if(find(labels.begin(),labels.end(),lb)!=labels.end()) {
        int count = 0;
        for(i=0; i<C.size(); ++i) {
          if(C[i]==lb)
            count++;
        }
        if(count>max_size) {
          return true;
        }
        labels.push_back(lb);
      }
    }
  }*/

  //return false;
}

