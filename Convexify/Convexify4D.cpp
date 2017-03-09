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

vector<float>
takeAverage(vector<CParticle4D>& P)
{
	vector<float> av(4, 0);
	for (int i = 0; i < P.size(); ++i)
	{
		av[0] += P[i].m_X;
		av[1] += P[i].m_Y;
		av[2] += P[i].m_Z;
		av[3] += P[i].m_T;
	}
	for (int j = 0; j < av.size(); ++j)
	{
		av[j] /= (float)P.size();
	}
	return av;
}

float norm(vector<float>& v)
{
	float sum = 0;
	for (int i = 0; i < v.size(); ++i)
	{
		sum += v[i] * v[i];
	}
	return sqrt(sum);
}

float det3(float a, float b, float c, float d, float e, float f, float g, float h, float i)
{
	return a* e* i - a* f* h - b* d* i + b* f* g + c* d* h - c* e* g;
}

vector<float>
normalVector(CParticle4D points[])
{
	float eps = 1.0e-6;
	float x1 = points[0].m_X, x2 = points[1].m_X, x3 = points[2].m_X, x4 = points[3].m_X;
	float y1 = points[0].m_Y, y2 = points[1].m_Y, y3 = points[2].m_Y, y4 = points[3].m_Y;
	float z1 = points[0].m_Z, z2 = points[1].m_Z, z3 = points[2].m_Z, z4 = points[3].m_Z;
	float t1 = points[0].m_T, t2 = points[1].m_T, t3 = points[2].m_T, t4 = points[3].m_T;
	/*float a = det3(y2 - y1, z2 - z1, t2 - t1, y3 - y2, z3 - z2, t3 - t2, y4 - y3, z4 - z3, t4 - t3);
	float b = -det3(x2 - x1, z2 - z1, t2 - t1, x3 - x2, z3 - z2, t3 - t2, x4 - x3, z4 - z3, t4 - t3);
	float c = det3(x2 - x1, y2 - y1, t2 - t1, x3 - x2, y3 - y2, t3 - t2, x4 - x3, y4 - y3, t4 - t3);
	float d = -det3(x2 - x1, y2 - y1, z2 - z1, x3 - x2, y3 - y2, z3 - z2, x4 - x3, y4 - y3, z4 - z3);*/
	float a = det3(y2 - y1, z2 - z1, t2 - t1, y3 - y1, z3 - z1, t3 - t1, y4 - y1, z4 - z1, t4 - t1);
	float b = -det3(x2 - x1, z2 - z1, t2 - t1, x3 - x1, z3 - z1, t3 - t1, x4 - x1, z4 - z1, t4 - t1);
	float c = det3(x2 - x1, y2 - y1, t2 - t1, x3 - x1, y3 - y1, t3 - t1, x4 - x1, y4 - y1, t4 - t1);
	float d = -det3(x2 - x1, y2 - y1, z2 - z1, x3 - x1, y3 - y1, z3 - z1, x4 - x1, y4 - y1, z4 - z1);

	float len = sqrt(a*a + b*b + c*c + d*d);
	vector<float> vec(4, 0);
	if (len > 0)
	{
		vec[0] = a / len; vec[1] = b / len; vec[2] = c / len; vec[3] = d / len;
	}
	return vec;
}

bool
positiveSide(vector<float>& outer, float x, float y, float z, float t)
{
	return outer[0] * x + outer[1] * y + outer[2] * z + outer[3] * t > 0;
}


vector<float>
outerVector(vector<float>& normal, vector<float>& center, CParticle4D& p0)
{
	if (positiveSide(normal, center[0] - p0.m_X, center[1] - p0.m_Y, center[2] - p0.m_Z, center[3] - p0.m_T))
	{
		for (int i = 0; i < normal.size(); ++i)
		{
			normal[i] = -normal[i];
		}
	}
	return normal;
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
	CParticle4D top_left;
	CParticle4D bottom_right;
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
			CParticle4D p4(x, y, z, t);
			if (i == 0)
			{
				top_left = p4;
				bottom_right = p4;
			}
			else
			{
				top_left.m_X = Min(top_left.m_X, p4.m_X);
				top_left.m_Y = Min(top_left.m_Y, p4.m_Y);
				top_left.m_Z = Min(top_left.m_Z, p4.m_Z);
				top_left.m_T = Min(top_left.m_T, p4.m_T);
				bottom_right.m_X = Max(bottom_right.m_X, p4.m_X);
				bottom_right.m_Y = Max(bottom_right.m_Y, p4.m_Y);
				bottom_right.m_Z = Max(bottom_right.m_Z, p4.m_Z);
				bottom_right.m_T = Max(bottom_right.m_T, p4.m_T);
			}
			P.push_back(p4);
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
	//load figure-ground segmentation
	int ndim;
	const int* dims;
	vector<unsigned char> L;
	{
		mxClassID classL;
		LoadData(L, prhs[2], classL, ndim, &dims);
	}
	assert(ndim == 4);
	//make sure the bounding box is within the target foreground.
	top_left.m_X = Max(top_left.m_X, 0);
	top_left.m_Y = Max(top_left.m_Y, 0);
	top_left.m_Z = Max(top_left.m_Z, 0);
	top_left.m_T = Max(top_left.m_T, 0);
	bottom_right.m_X = Min(bottom_right.m_X, dims[0] - 1);
	bottom_right.m_Y = Min(bottom_right.m_Y, dims[1] - 1);
	bottom_right.m_Z = Min(bottom_right.m_Z, dims[2] - 1);
	bottom_right.m_T = Min(bottom_right.m_T, dims[3] - 1);

	int nvoxels = numberOfElements(ndim, dims);
	printf("ndim = %d, dims=[%d, %d, %d, %d], nvoxels=%d, %d\n",
		ndim, dims[0], dims[1], dims[2], dims[3], nvoxels, L.size());
	printf("Bounding box =(%d, %d, %d, %d) - (%d, %d, %d, %d)\n",
		top_left.m_X, top_left.m_Y, top_left.m_Z, top_left.m_T, bottom_right.m_X, bottom_right.m_Y, bottom_right.m_Z, bottom_right.m_T);

	//clear voxels outside the bounding box
	for (int m = 0; m < dims[3]; ++m)
	{
		for (int i = 0; i < dims[2]; ++i)
		{
			for (int j = 0; j < dims[1]; ++j)
			{
				for (int k = 0; k < dims[0]; ++k)
				{
					if (m<top_left.m_T || m>bottom_right.m_T || i<top_left.m_Z || i>bottom_right.m_Z ||
						j<top_left.m_Y || j>bottom_right.m_Y || k<top_left.m_X || k>bottom_right.m_X)
					{
						SetData4(L, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0);
					}
				}
			}
		}
	}

	vector<float> center = takeAverage(P);
	for (int n = 0; n < dimsI[0]; ++n)
	{
		//printf("%d\n", n);
		CParticle4D points[4];
		for (int i = 0; i < 4; ++i)
		{
			int idx = GetData2(I, n, i, dimsI[0], dimsI[1], -1);
			points[3-i] = P[idx];
		}
		vector<float> normal = normalVector(points);
		vector<float> outer = outerVector(normal, center, points[0]);
		/*printf("(%d,%d,%d,%d),(%d,%d,%d,%d),(%d,%d,%d,%d),(%d,%d,%d,%d)=>(%2.2f,%2.2f,%2.2f,%2.2f)\n",
			(points[0].m_X - 5) / 10, (points[0].m_Y - 5) / 10, (points[0].m_Z - 5) / 10, (points[0].m_T - 5) / 10,
			(points[1].m_X - 5) / 10, (points[1].m_Y - 5) / 10, (points[1].m_Z - 5) / 10, (points[1].m_T - 5) / 10,
			(points[2].m_X - 5) / 10, (points[2].m_Y - 5) / 10, (points[2].m_Z - 5) / 10, (points[2].m_T - 5) / 10,
			(points[3].m_X - 5) / 10, (points[3].m_Y - 5) / 10, (points[3].m_Z - 5) / 10, (points[3].m_T - 5) / 10,
			outer[0], outer[1], outer[2], outer[3]);*/
		if (norm(outer) <= 0)
		{
			//printf("Zero outer normal vector.");
			continue;
		}
		for (int m = top_left.m_T; m <= bottom_right.m_T; ++m)
		{
			float t2 = m - points[0].m_T;
			for (int i = top_left.m_Z; i <= bottom_right.m_Z; ++i)
			{
				float z2 = i - points[0].m_Z;
				for (int j = top_left.m_Y; j<= bottom_right.m_Y; ++j)
				{
					float y2 = j - points[0].m_Y;
					for (int k = top_left.m_X; k<= bottom_right.m_X; ++k)
					{
						float x2 = k - points[0].m_X;
						if (GetData4(L, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)1))
						{
							if (positiveSide(outer, x2, y2, z2, t2)) //outside convex hull
							{
								SetData4(L, k, j, i, m, dims[0], dims[1], dims[2], dims[3], (unsigned char)0);
							}
						}
					}
				}
			}
		}
	}

	//printf("Done\n");
	if (nlhs >= 1)
	{
		plhs[0] = StoreData(L,mxUINT8_CLASS,ndim,dims);
	}
	mexUnlock();
}

