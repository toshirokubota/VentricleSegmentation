#ifndef ___SZ_MEDIAL_AXIS_H___
#define ___SZ_MEDIAL_AXIS_H___
#include <szParticle.h>
#include <vector>
using namespace std;

struct MedialAxis
{
	MedialAxis(int x=0, int y=0, int z=0): Contact(x, y, z) {}
	MedialAxis(const CParticle& p): Contact(p){}

	bool operator == (const MedialAxis& p) const
	{
		if(Contact == p.Contact)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	vector<CParticle> Core;
	CParticle Contact;
};

#endif /* ___SZ_MEDIAL_AXIS_H___ */