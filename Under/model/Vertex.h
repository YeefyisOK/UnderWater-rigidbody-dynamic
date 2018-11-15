#pragma once

#include "Point3D.h"

class CMesh;

class CVertex
	:public CPoint3D
{
public:
	CVertex();
	~CVertex();

	CVertex &operator= (const CPoint3D &pt);

	void SetId(int id)
	{
		m_id = id;
	}

	int GetId()
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

	void SetTexCoord(double uv[2])
	{
		m_uv[0] = uv[0];
		m_uv[1] = uv[1];
	}
	double * GetTexCoord()
	{
		return m_uv;
	}


	void SetNormal(CVector3D &normal)
	{
		m_normal = normal;
		m_normal.Normalize();
	}
	CVector3D &GetNormal()
	{
		return m_normal;
	}

	void AddNearEdge(int eId)
	{
		m_vecNearEdges.push_back(eId);
	}

	std::vector<int> &GetNearEdges()
	{
		return m_vecNearEdges;
	}

	bool GetAdjacentFacesID(std::vector<int>& adjFaces) const;
	bool GetAdjacentVerticesId(std::vector<int> &adjVertices) const;

	const std::vector<int> & GetAdjacentEdgesId() const
	{
		return m_vecNearEdges;
	}
private:
	CMesh *m_pParent;

	CVector3D m_normal;

	std::vector<int> m_vecNearEdges;

	int m_id;
	double m_uv[2];
};

