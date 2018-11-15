#include "stdafx.h"
#include "face.h"
#include "Mesh.h"


CFace::CFace()
{
	m_pParent = nullptr;
	m_id = -1;
}


CFace::~CFace()
{
}

void CFace::ComputeNormal()
{
	CPoint3D &v0 = m_pParent->GetVertexAt(m_vId[0]);
	CPoint3D &v1 = m_pParent->GetVertexAt(m_vId[1]);
	CPoint3D &v2 = m_pParent->GetVertexAt(m_vId[2]);

	m_normal = (v1 - v0) * (v2 - v0);
	m_normal.Normalize();
}

bool CFace::GetAdjacentFacesId(int fId[3]) const
{
	for (int i = 0; i < 3; ++i)
	{
		CEdge &edge = m_pParent->GetEdgeAt(m_edgeId[i]);
		fId[i] = edge.GetAnotherFaceId(m_id);
	}
	return true;
}

bool CFace::GetAdjacentFacesId(std::vector<int> &adjFaces) const
{
	adjFaces.clear();

	std::set<int> adjFacesSet;
	for (int i = 0; i < 3; ++i)
	{
		CVertex &vertex = m_pParent->GetVertexAt(m_vId[i]);
		const std::vector<int> &adjEdgesId = vertex.GetAdjacentEdgesId();

		for (size_t j = 0; j < adjEdgesId.size(); ++j)
		{
			CEdge &edge = m_pParent->GetEdgeAt(adjEdgesId[j]);
			int fId = edge.GetAnotherFaceId(m_id);

			if (fId >= 0)
			{
				adjFacesSet.insert(fId);
			}
		}
	}

	//convert to vector
	std::copy(adjFacesSet.begin(), adjFacesSet.end(), std::back_inserter(adjFaces));
	return true;
}