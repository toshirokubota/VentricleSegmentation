#include <szExtractNodule.h>
#include <szDefaultParam.h>
#include <szMexUtility.h>
#include <szIsFunctions.h>
#include <szCutOutNodule.h>
#include <szConvexHull3D.h>
#include <szMiscOperations.h>
#include <szMyNeighborOp.h>

int
doNoduleExtraction(const vector<unsigned char>& L, //input - foreground
				   const vector<float>& D, //input - distance map
				   //IMPORTANT! - B and B0 need to be initialized to zero by the calling function
				   vector<unsigned char>& B, //output - nodule segmentation 
				   vector<unsigned char>& B0, //output - nodule segmentation after the 1st phase
				   const vector<int>& vseeds, //input - seed positions in a set of triplets (x, y, z)
				   bool bConvexHull, //enable/disable a convex hull operation
				   bool bRefineSeg, //enable/disable an operation that refines segmentation boundary
				   int ndim,
				   const int* dims)
{
	int nvoxels = numberOfElements(ndim,dims);
	int attachType = SZ_LABEL_NO_NODULE;
	if(IsSolitary(L, vseeds[0], vseeds[1], vseeds[2], nvoxels/8, dims))
	{
		printf("Solitary nodule.\n");
		B = L;
		selectCluster(B, vseeds[0], vseeds[1], vseeds[2], dims);
		attachType = SZ_LABEL_SOLITARY_NODULE;
	}
	else
	{
		attachType = SZ_LABEL_JUXTAPOSITION_NODULE;

		printf("Non-Solitary nodule.\n");
		for(int i=0; i<vseeds.size(); i+=3)
		{
			vector<unsigned char> B2(nvoxels,0);
			CParticle seed(vseeds[i], vseeds[i+1], vseeds[i+2], 0);
			vector<CParticle> vseeds(1, seed);
			vseeds = TraceNew2(B2, L, D, vseeds, 0.8, NonStrictDownHillTrace, dims);
			MergeLabels(B0, B2, dims);
			vseeds = TraceNew2(B2, L, D, vseeds, 1.01, StrictDownHillTrace, dims);
			MergeLabels(B, B2, dims);
		}
		selectCluster(B, vseeds[0], vseeds[1], vseeds[2], dims);
		if(bConvexHull)
		{
			vector<unsigned char> Surf(nvoxels,0);
			vector<unsigned char> CH(nvoxels,0);
			vector<unsigned char> Sg(nvoxels,0);
			preConvexHull(Surf,B,dims);
			bool bRet = ConvexHull3D (CH,Surf,ndim,dims); 
			if(!bRet)
				printf("Convex Hull failed...\n");
			else
				postConvexHull(B, CH, L, dims);
		}
		if(bRefineSeg)
		{
			vector<int> causalNbhd = MakeCausalNeighborhood(3, dims);
			B = refineSegmentation(B, L, causalNbhd, 3, dims);
		}
	}
	return attachType;
}