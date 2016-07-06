#include <szParticle.h>

template <class T>
bool
growLE(const CParticle& p,
	 const vector<T>& D, 
	 const int* dims)
{
	float dval = GetData3(D, p.m_X, p.m_Y, p.m_Z, dims[0], dims[1], dims[2], (T)0);
	vector<CParticle> region(1, p);
	vector<CParticle> vp(1, p);
	while(vp.empty() == false)
	{
		vector<CParticle> vp2;
		for(int i=0; i<vp.size(); ++i)
		{
			int x = vp[i].m_X;
			int y = vp[i].m_Y;
			int z = vp[i].m_Z;
			for(int m=0; m<NumNeighbors; ++m)
			{
				int x2 = x + XOffset[m];
				int y2 = y + YOffset[m];
				int z2 = z + ZOffset[m];
				float dval2 = GetData3(D, x2, y2, z2, dims[0], dims[1], dims[2], (T)0);
				if(dval2 > dval)
				{
					return false;
				}
				else if(dval2 == dval)
				{
					CParticle q(x2, y2, z2);
					if(find(region.begin(), region.end(), q) == region.end())
					{
						vp2.push_back(q);
						region.push_back(q);
					}
				}
			}
		}
		vp = vp2;
	}
	return true;
}

/*
Spurious local maximum is a local maximum component that is adjacent to a component with a 
higher distance value.  It was picked as a local maximum due to non-strict condition.
*/
template <class T>
void
RemoveSpurious(vector<unsigned char>& P,
		 const vector<T>& D,
		 const int* dims)
{
	vector<CParticle> vp;
	for(int i=0; i<dims[2]; ++i)
	{
		for(int j=0; j<dims[1]; ++j)
		{
			for(int k=0; k<dims[0]; ++k)
			{
				if(GetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
				{
					CParticle p(k, j, i);
					if(growLE(p, D, dims) == false)
					{
						SetData3(P, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0);
					}
				}
			}
		}
	}
}

/* 
*/
template<class Item>
void
LocalMaximum(vector<unsigned char>& M, 
             const vector<Item>& V, 
             const vector<unsigned char>& L, //ROI mask
             const vector<int>& nbh,
             bool strict,
             int ndim,
             const int* dims) {

  int nvoxels = numberOfElements(ndim,dims);
  if(M.size()<nvoxels || V.size()<nvoxels || L.size()<nvoxels) {
    //mexPrintf("An input image does not contain enough data.\n");
    return;
  }

  int i;
  //construct a subscript representation of the neighborhood
  //for efficient boundary check
  vector<int> vcoords;
  for(i=0; i<nbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
    vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
  }

  for(i=0; i<nvoxels; ++i) {
    unsigned char lb=L[i];
    if(!lb)
	{
      SetData(M,i,(unsigned char) 0);
      continue; //not inside ROI
	}
    bool bLM=true;
    Item origV=V[i];
    
    vector<int> vsub = Ind2Sub(i,ndim,dims);
    for(int j=0; j<nbh.size(); ++j) {
      if(NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims)) { 
        int k=i+nbh[j];
        if(L[k]) {
          Item v2=V[k];
          if(strict) {// strictly maximum
            if(v2>=origV) {
              bLM=false;
              break;
            }
          }
          else {
            if(v2>origV) {
              bLM=false;
              break;
            }
          }
        }
      }
    }
        
    if(bLM) {
      SetData(M,i,(unsigned char) 1);
    }
    else {
      SetData(M,i,(unsigned char) 0);
    }
  }
}


/* 
*/
template<class Item>
void
LocalMinimum(vector<unsigned char>& M, 
             const vector<Item>& V, 
             const vector<unsigned char>& L,
             const vector<int>& nbh,
             bool strict,
             int ndim,
             const int* dims) {
  int nvoxels = numberOfElements(ndim,dims);
  if(M.size()<nvoxels || V.size()<nvoxels || L.size()<nvoxels) {
    //mexPrintf("An input image does not contain enough data.\n");
    return;
  }

  int i;
  //construct a subscript representation of the neighborhood
  //for efficient boundary check
  vector<int> vcoords;
  for(i=0; i<nbh.size(); ++i) {
    vector<int> vsub = Ind2SubCentered(nbh[i],ndim,dims);
    vcoords.insert(vcoords.end(),vsub.begin(),vsub.end());
  }

  for(i=0; i<nvoxels; ++i) {
    unsigned char lb=L[i];
    if(!lb)
      continue; //not inside ROI
    bool bLM=true;
    bool b0;
    Item origV=V[i];
    
    vector<int> vsub = Ind2Sub(i,ndim,dims);
    for(int j=0; j<nbh.size(); ++j) {
      if(NeighborCheck(vsub.begin(),vcoords.begin()+j*ndim,ndim,dims)) { 
        int k=i+nbh[j];
        if(L[k]) {
          bool b;
          Item v2=V[k];
          if(strict) {// strictly maximum
            if(v2<=origV) {
              bLM=false;
              break;
            }
          }
          else {
            if(v2<origV) {
              bLM=false;
              break;
            }
          }
        }
      }
    }
        
    if(bLM) {
      SetData(M,i,(unsigned char) 1);
    }
    else {
      SetData(M,i,(unsigned char) 0);
    }
  }
}

