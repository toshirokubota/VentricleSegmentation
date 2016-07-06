#ifndef _SZ_DMAP_SEGMENTATION_H_
#define _SZ_DMAP_SEGMENTATION_H_

#include<vector>
using namespace std;

void
initialTrace(vector<unsigned char>& S,
          vector<float>& D,
          vector<unsigned char>& M,
          int x, int y, int z,
          float lowThres,
          float highThres,
          vector<int>& seeds,
          const int* dims);

bool
initialTraceWall(vector<unsigned char>& S,
                 vector<float>& D,
                 vector<unsigned char>& M, //intermediate trace
                 int x, int y, int z,
                 float lowThres,
                 float highThres,
                 vector<unsigned char>& wall,
                 vector<int>& seeds,
                 const int* dims);

void
DistanceMapSegmentation(vector<unsigned char>& S,   //OUTPUT - segmentation
                        vector<int>& peaks,         //seeds for the trace
                        vector<float>& D,           //INPUT - distance map
                        vector<unsigned char>& L,   //INPUT - foreground segmentation
                        float thres,                //INPUT - fitness threshold
                        const int* dims);

#endif /* _SZ_DMAP_SEGMENTATION_H_ */