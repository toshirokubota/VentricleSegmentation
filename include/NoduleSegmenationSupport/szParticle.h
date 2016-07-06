#ifndef ___SZ_CPARTICLE_H___
#define ___SZ_CPARTICLE_H___

class CParticle
{
public:
	CParticle(int x=0, int y=0, int z=0, float life=0, int x0=-1, int y0=-1, int z0=-1)
	{
		m_X = x;
		m_Y = y;
		m_Z = z;
		if(x0<0)
		{
			m_X0 = x;
		}
		else
		{
			m_X0 = x0;
		}
		if(y0<0)
		{
			m_Y0 = y;
		}
		else
		{
			m_Y0 = y0;
		}
		if(z0<0)
		{
			m_Z0 = z;
		}
		else
		{
			m_Z0 = z0;
		}
		m_Life = life;
	}
	bool operator == (const CParticle& p) const
	{
		if(m_X == p.m_X && m_Y == p.m_Y && m_Z == p.m_Z)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator < (const CParticle& p) const
	{
		if(this->m_Life < p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator <= (const CParticle& p) const
	{
		if(this->m_Life <= p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator > (const CParticle& p) const
	{
		if(this->m_Life > p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
	bool operator >= (const CParticle& p) const
	{
		if(this->m_Life >= p.m_Life)
		{
			return true;
		}
		else
		{
			return false;
		}
	}
public:
	int m_X;
	int m_Y;
	int m_Z;
	int m_X0;
	int m_Y0;
	int m_Z0;
	float m_Life;
};

#endif /* ___SZ_CPARTICLE_H___ */