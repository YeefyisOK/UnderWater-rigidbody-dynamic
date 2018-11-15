#pragma once
#include "Point3D.h"
class CBBox
{
public:
	CBBox( );
	CBBox(const CPoint3D &ptMin, const CPoint3D &ptMax)
	{
		m_ptMin = ptMin;
		m_ptMax = ptMax;
	}
	~CBBox();

	double DiagonalLength() const
	{
		return m_ptMin.DistanceWith(m_ptMax);
	}

	double SquareDiagonalLength() const
	{
		return m_ptMin.SquareDistanceWith(m_ptMax);
	}

	void Set(const CPoint3D &ptMin, const CPoint3D &ptMax)
	{
		m_ptMin = ptMin;
		m_ptMax = ptMax;
	}

	const CPoint3D &MinPoint() const
	{
		return m_ptMin;
	}

	const CPoint3D &MaxPoint() const
	{
		return m_ptMax;
	}

	void SetEmpty()
	{
		m_ptMin.Set(0, 0, 0);
		m_ptMax.Set(0, 0, 0);
	}

	void MergeWith(const CBBox &bbox);

	void Expan(const CBBox &bbox);

	void Translate(const CVector3D &offset);

	CPoint3D Center() const
	{
		return (m_ptMin + m_ptMax) * 0.5;
	}

	static CBBox GetBBox(std::vector<CPoint3D> &vecPnts)
	{
		if (vecPnts.empty())
		{
			return CBBox();
		}
		return GetBBox(&(vecPnts.front()), vecPnts.size());
	}
	static CBBox GetBBox(const CPoint3D *pPnts, size_t nPnts);
private:
	CPoint3D m_ptMin, m_ptMax;
};

