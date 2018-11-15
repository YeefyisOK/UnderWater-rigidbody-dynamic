#pragma once
#include "Primitive.h"
#include "Point3D.h"

class CMesh;

class CFace
	:public CPrimitive
{
public:
	CFace();
	~CFace();

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

	void SetVerticeId(int v0, int v1, int v2)
	{
		m_vId[0] = v0;
		m_vId[1] = v1;
		m_vId[2] = v2;
	}
	void SetVerticeId(int vId[3])
	{
		SetVerticeId(vId[0], vId[1], vId[2]);
	}

	int* GetVerticeId()
	{
		return m_vId;
	}

	void SetEdgesId(int e0, int e1, int e2)
	{
		m_edgeId[0] = e0;
		m_edgeId[1] = e1;
		m_edgeId[2] = e2;
	}
	void SetEdgesId(int eId[3])
	{
		SetEdgesId(eId[0], eId[1], eId[2]);
	}
	int *GetEdgesId()
	{
		return m_edgeId;
	}

	void ComputeNormal();

	CVector3D &GetNormal()
	{
		return m_normal;
	}

	bool GetAdjacentFacesId(int fId[3]) const;
	bool GetAdjacentFacesId(std::vector<int> &adjFaces) const;
private:
	int m_vId[3];
	int m_edgeId[3];
	CVector3D m_normal;

	int m_id;

	CMesh *m_pParent;
};

