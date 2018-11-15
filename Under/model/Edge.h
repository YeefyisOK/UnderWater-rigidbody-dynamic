#pragma once

#include "Primitive.h"

class CMesh;

class CEdge
	:public CPrimitive
{
public:
	CEdge(int v0 = -1, int v1 = -1);
	~CEdge();

	void SetId(int id)
	{
		m_id = id;
	}

	int GetId() const
	{
		return m_id;
	}

	void SetParent(CMesh *pParent)
	{
		m_pParent = pParent;
	}

	CMesh *GetParent()
	{
		return m_pParent;
	}

	int GetAnotherVertexId(int vId) const
	{
		if (m_vId[0] == vId)
		{
			return m_vId[1];
		}
		return m_vId[0];
	}

	int GetAnotherFaceId(int fId) const
	{
		if (m_left == fId)
		{
			return m_right;
		}
		return m_left;
	}

	void AddFace(int fId)
	{
		if (m_left == -1)
		{
			m_left = fId;
		} 
		else
		{
			m_right = fId;
		}
	}

	void GetAdjacentFaceId(int &left, int &right) const
	{
		left = m_left;
		right = m_right;
	}

	int *GetVerticeId()
	{
		return m_vId;
	}
private:
	int m_vId[2];

	int m_left, m_right;

	int m_id;
	CMesh *m_pParent;
};

