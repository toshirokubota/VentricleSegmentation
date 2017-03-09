#include <iostream>
using namespace std;
#include <stdio.h>

#ifdef MEX_DLL
#include <mex.h>
#include "mexFileIO.h"
#endif

#include <szMexUtility.h>
#include <szMexUtilityTemplate.h>
#include <szDefaultParam.h>
#include <szConvexHull3D.h>
#include <szCutOutNodule.h>
#include <szMiscOperations.h>
#include <szParticle4D.h>

float det2(float a, float b, float c, float d)
{
	return a * d - c * b;
}

vector<float>
outerVector(CParticle4D points[])
{
	float eps = 1.0e-6;
	float x1 = points[0].m_X, x2 = points[1].m_X, x3 = points[2].m_X;
	float y1 = points[0].m_Y, y2 = points[1].m_Y, y3 = points[2].m_Y;
	float z1 = points[0].m_Z, z2 = points[1].m_Z, z3 = points[2].m_Z;
	float a = -det2(y2 - y1, z2 - z1, y3 - y2, z3 - z2);
	float b = det2(x2 - x1, z2 - z1, x3 - x2, z3 - z2);
	float c = -det2(x2 - x1, y2 - y1, x3 - x2, y3 - y2);

	float len = sqrt(a*a + b*b + c*c);
	vector<float> vec(3, 0);
	if (len > 0)
	{
		vec[0] = a / len; vec[1] = b / len; vec[2] = c / len;
	}
	return vec;
}


bool
positiveSide(vector<float>& outer, float x, float y, float z)
{
	return outer[0] * x + outer[1] * y + outer[2] * z > 0;
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	if (nrhs < 1 || nlhs < 1)
	{
		mexErrMsgTxt("Usage: B = FigureGroundSeparation(A, [L])");
		return;
	}

	//load points
	vector<CParticle4D> P;
	{
		vector<int> P0;
		int ndimP;
		const int* dimsP;
		mxClassID classP;
		LoadData(P0, prhs[0], classP, ndimP, &dimsP);
		for (int i = 0; i < dimsP[0]; ++i)
		{
			int x = GetData2(P0, i, 0, dimsP[0], dimsP[1], 0);
			int y = GetData2(P0, i, 1, dimsP[0], dimsP[1], 0);
			int z = GetData2(P0, i, 2, dimsP[0], dimsP[1], 0);
			int t = GetData2(P0, i, 3, dimsP[0], dimsP[1], 0);
			P.push_back(CParticle4D(x, y, z, t));
		}
	}
	//load indices
	vector<int> I;
	const int* dimsI;
	{
		int ndimI;
		mxClassID classI;
		LoadData(I, prhs[1], classI, ndimI, &dimsI);
	}
	//load data size
	vector<int> sz;
	{
		int ndimS;
		const int* dimsS;
		mxClassID classS;
		LoadData(sz, prhs[2], classS, ndimS, &dimsS);
	}
	int ndim = 3;
	int dims[] = { sz[0], sz[1], sz[2]};
	int nvoxels = numberOfElements(ndim, dims);
	vector<unsigned char> B(nvoxels, 1);
	printf("ndim = %d, dims=[%d, %d, %d, %d], nvoxels=%d, %d\n", 
		ndim, dims[0], dims[1], dims[2], dims[3], nvoxels, B.size());
	for (int n = 0; n < dimsI[0]; ++n)
	{
		printf("%d\n", n);
		CParticle4D points[3];
		for (int i = 0; i < 3; ++i)
		{
			int idx = GetData2(I, n, i, dimsI[0], dimsI[1], -1);
			points[i] = P[idx];
		}
		vector<float> outer = outerVector(points);
		printf("(%d,%d,%d),(%d,%d,%d),(%d,%d,%d)=>(%2.2f,%2.2f,%2.2f)\n",
			(points[0].m_X - 5) / 10, (points[0].m_Y - 5) / 10, (points[0].m_Z - 5) / 10,
			(points[1].m_X - 5) / 10, (points[1].m_Y - 5) / 10, (points[1].m_Z - 5) / 10,
			(points[2].m_X - 5) / 10, (points[2].m_Y - 5) / 10, (points[2].m_Z - 5) / 10,
			outer[0], outer[1], outer[2]);
		for (int i = 0; i < dims[2]; ++i)
		{
			float z2 = i - points[0].m_Z;
			for (int j = 0; j < dims[1]; ++j)
			{
				float y2 = j - points[0].m_Y;
				for (int k = 0; k < dims[0]; ++k)
				{
					float x2 = k - points[0].m_X;
					if (GetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0))
					{
						if (positiveSide(outer, x2, y2, z2))
						{
							SetData3(B, k, j, i, dims[0], dims[1], dims[2], (unsigned char)0);
						}
					}
				}
			}
		}
	}

	printf("Done\n");
	if (nlhs >= 1)
	{
		plhs[0] = StoreData(B,mxUINT8_CLASS,ndim,dims);
	}
	mexUnlock();
}

