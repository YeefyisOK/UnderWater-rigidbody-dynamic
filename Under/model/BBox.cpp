#include "stdafx.h"
#include "BBox.h"


CBBox::CBBox()
{
}


CBBox::~CBBox()
{
}


void CBBox::MergeWith(const CBBox &bbox)
{
	m_ptMin[0] = std::min(m_ptMin[0], bbox.m_ptMin[0]);
	m_ptMin[1] = std::min(m_ptMin[1], bbox.m_ptMin[1]);
	m_ptMin[2] = std::min(m_ptMin[2], bbox.m_ptMin[2]);

	m_ptMax[0] = std::max(m_ptMax[0], bbox.m_ptMax[0]);
	m_ptMax[1] = std::max(m_ptMax[1], bbox.m_ptMax[1]);
	m_ptMax[2] = std::max(m_ptMax[2], bbox.m_ptMax[2]);
}

void CBBox::Translate(const CVector3D &offset)
{
	m_ptMin += offset;
	m_ptMax += offset;
}

void CBBox::Expan(const CBBox &bbox)
{
	for (int i = 0; i < 3; ++i)
	{
		if (bbox.m_ptMin[i] < 0)
		{
			m_ptMin[i] += bbox.m_ptMin[i];
		}

		if (bbox.m_ptMax[i] > 0)
		{
			m_ptMax[i] += bbox.m_ptMax[i];
		}
	}
}

CBBox CBBox::GetBBox(const CPoint3D *pPoints, size_t nPnts)
{
	CBBox bbox;
	if (nPnts == 0)
	{
		return bbox;
	}

	bbox.m_ptMin = bbox.m_ptMax = pPoints[0];
	for (size_t i = 1; i < nPnts; ++i)
	{
		//x
		if (pPoints[i][0] < bbox.m_ptMin[0])
		{
			bbox.m_ptMin[0] = pPoints[i][0];
		} 
		else if (pPoints[i][0] > bbox.m_ptMax[0])
		{
			bbox.m_ptMax[0] = pPoints[i][0];
		}

		//y
		if (pPoints[i][1] < bbox.m_ptMin[1])
		{
			bbox.m_ptMin[1] = pPoints[i][1];
		}
		else if (pPoints[i][1] > bbox.m_ptMax[1])
		{
			bbox.m_ptMax[1] = pPoints[i][1];
		}

		//z
		if (pPoints[i][2] < bbox.m_ptMin[2])
		{
			bbox.m_ptMin[2] = pPoints[i][2];
		}
		else if (pPoints[i][2] > bbox.m_ptMax[2])
		{
			bbox.m_ptMax[2] = pPoints[i][2];
		}
	}
	return bbox;
}