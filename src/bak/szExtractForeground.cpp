
#ifdef MEX_DLL
#include <mex.h>
#endif
#include <szDefaultParam.h>
#include <szExtractForeground.h>
#include <szFigureGroundSeparation.h>
#include <szFindSeeds.h>
#include <szMiscOperations.h>
#include <szMyNeighborOp.h>
#include <szConnectedComponent.h>

//do figure-ground segmentation and obtain seed locations
int
doForegroundExtraction(const vector<unsigned short>& A, //input image
					   vector<unsigned char>& L, //output - foreground segmentation 
					   vector<float>& Q, //output - probability image
					   vector<float>& D, //output - distance map image
					   vector<float>& Sp, //output - sphericity image
					   vector<unsigned char>& P, //output - sphericity local maxima
					   vector<int>& vseeds, //output - seed points
					   vector<float>& vSph, //sphericity values at the seed points
					   const vector<int>& vnoduleTypes, //input - nodule types of interest
					   //what to follow are various parameters to the RD segmentation
					   int RDIter,
					   float RDlambda,
					   //what to follow are various parameters to the seed search
					   int findSeedWidth, 
					   int findSeedCoreWidth,
					   float findSeedThres,
					   //a parameter that enables/disables persistent segmentation
					   bool bPersistent,
					   //dimension of the data
					   int ndim,
					   const int* dims)
{
	int noduleType = SZ_LABEL_NO_NODULE;
	for(int ti=0; ti<vnoduleTypes.size(); ++ti)
	{ 
		vector<unsigned char> Ls(A.size(),(unsigned char)0);
		vector<float> Qs(A.size(),(float)0);
		vector<float> Ds(A.size(),(float)0);
		vector<float> Sps(A.size(),(float)0);
		vector<unsigned char> Ps(A.size(),(unsigned char)0);
		vector<int> vSeedsQ;
		vector<float> vSphQ;

		printf("doForegroundExtraction: running with mode %d.  Current best mode is %d\n", 
			vnoduleTypes[ti], noduleType);

		FigureGroundSeparation(
			A,
			Ls,
			Qs,
			vnoduleTypes[ti],
			RDIter,
			RDlambda,
			dims);

		FillBackgroundHoles(Ls, 15, 1, dims);

		vSeedsQ = FindSeedsThresNew(Ls, Ds, Sps, Ps, vSphQ,
			findSeedWidth, 
			findSeedCoreWidth,
			findSeedThres,
			1,
			dims);
		if(vSeedsQ.empty() && bPersistent)
		{ 
			int model = -1;
			if(vnoduleTypes[ti] == SZ_LABEL_SOLID_NODULE)
			{
				printf("Trying with SZ_LABEL_SOLID_NODULE_LOW...\n");
				model = SZ_LABEL_SOLID_NODULE_LOW;
			}
			else if(vnoduleTypes[ti] == SZ_LABEL_NONSOLID_NODULE)
			{
				printf("Trying with SZ_LABEL_NONSOLID_NODULE_LOW...\n");
				model = SZ_LABEL_NONSOLID_NODULE_LOW;
			}
			if(model>0) 
			{				
				FigureGroundSeparation(
					A,
					Ls,
					Qs,
					model,
					RDIter,
					RDlambda,
					dims);	  
				vSeedsQ = FindSeedsThres(Ls, Ds, Sps, Ps, vSphQ,
					findSeedWidth, 
					findSeedCoreWidth,
					findSeedThres,
					dims);
			}
		}
		if(vnoduleTypes.size() == 1 ||
			compareForegroundFitness2(vSph, noduleType, vSphQ, vnoduleTypes[ti]) < 0)
		{
			L = Ls;
			D = Ds;
			P = Ps;
			Sp = Sps;
			Q = Qs;
			vseeds = vSeedsQ;
			vSph = vSphQ;
			noduleType = vnoduleTypes[ti];
		}
	}

	return noduleType;
}

