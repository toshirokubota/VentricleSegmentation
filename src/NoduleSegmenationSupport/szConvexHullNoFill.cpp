#ifdef MEX_DLL
#include <mex.h>
#endif
#include <fstream>
#include <vector>
#include <list>
#include <algorithm>
using namespace std;
#include <stdio.h>
#include <szDefaultParam.h>
#include <szMiscOperations.h>
#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szGeometry.h>
#include <szConvexHulFillNew.h>

/*
This is the meat of the ConvexHull routine.  
The input to the routine is a vector of vertices (CParticle) and
the result is retured also as a vector of vertices (CParticle).
The implementation is applicable only to 3D volume data and is based on the book:
M. de Berg, M. van Kreveld, M. Overmars, and O. Schwarzkopf,
"Computational Geometry - Algorithms and Applications", Springer.
The algorithm is described in page 236-239.

OUTPUT:
: result as a half-plane segmentation
INPUT:
vp: a set of vertices
*/

vector<CParticle>
ConvexHull3D (vector<CParticle>& vp)	 
{
	vector<CParticle> vpResult;
	if(vp.size() < 4)
	{
		printf("ConvexHull3D: At least 4 positions are required.  Only %d are found.\n", vp.size());
		return vpResult;  //return an empty vector
	}

	vector<szVertex*> v_pVertices; 
	vector<szFace*> v_pFaces;
	szConflictGraph CG;

	int nump = 0;
	int i;
	for(i=0; i<vp.size(); ++i) 
	{
		szVertex* pv = new szVertex;
		pv->setX(vp[i].m_X);
		pv->setY(vp[i].m_Y);
		pv->setZ(vp[i].m_Z);
		v_pVertices.push_back(pv);
		nump++;
	}
	//printf("ConvexHull3D: There are %d positions in the volume.\n", nump);

	RandomizeVertices(v_pVertices);

	bool bFlag = FindInitialHull(v_pVertices,v_pFaces);
	if(bFlag == false) 
	{
		printf("No non-coplaner vertices are found in %d voxels ...  No further processing is necessary for ConvexHull3D.\n", 
			v_pVertices.size());
		return false;
	}
	int numVerticesInHull = 4;
	InitializeConflictGraph(CG,v_pVertices,v_pFaces);

	bool bSuccess = true;
	bool track=false;
	for(i=numVerticesInHull; i<nump && bSuccess; ++i) 
	{
		szVertex* pv=*(v_pVertices.begin()+i);
		vector<szFace*> v_pConflictFaces=FindConflictingFaces(CG,pv);
		if(!v_pConflictFaces.empty()) 
		{//there are conflicts
			vector<szEdge*> v_pHorizontalEdges=FindHorizontalEdges(pv,v_pConflictFaces);
			vector<szFace*> v_pNewFaces = MakeNewVisibleFaces(pv,v_pHorizontalEdges);
			bSuccess=UpdateGeometricalPrimitives(pv,v_pNewFaces,v_pConflictFaces,v_pVertices,v_pFaces,track);
			UpdateConflictGraph(CG,pv,v_pNewFaces,v_pConflictFaces);
		}
		else 
		{
			//the point is inside the convex hull
		}
		//if(bSuccess)
		//  bSuccess=IntegrityCheck(v_pFaces);
	}

	if(!bSuccess || v_pFaces.size()<4) 
	{
		printf("Convex Hull failed. %d %d\n", bSuccess, v_pFaces.size());
	}
	else 
	{
		for(int i=0; i<v_pFaces.size(); ++i)
		{
			szVertex* pv = v_pFaces[n]->getOuterComponent()->getOrigin();
			int x = pv->getX();
			int y = pv->getY();
			int z = pv->getZ();
			vpResult.push_back(CParticle(x, y, z));
		}
	}

	for(i=0; i<v_pFaces.size(); ++i)
		delete v_pFaces[i];
	for(i=0; i<v_pVertices.size(); ++i)
		delete v_pVertices[i];

	return vpResult;  //return an empty vector
}

