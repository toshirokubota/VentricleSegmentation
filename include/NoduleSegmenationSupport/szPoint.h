#ifndef ___SZ_POINT_H___
#define ___SZ_POINT_H___

struct Point3D
{
	Point3D(int x=0, int y=0, int z=0)
	{
		X=x;
		Y=y;
		Z=z;
	}

	int X;
	int Y;
	int Z;
};

#endif /* ___SZ_POINT_H___ */