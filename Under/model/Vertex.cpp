#include "stdafx.h"
#include "Vertex.h"
#include "Edge.h"
#include "Mesh.h"

CVertex::CVertex()
{
	m_pParent = nullptr;
	m_id = -1;
}


CVertex::~CVertex()
{
}

CVertex &CVertex::operator=(const CPoint3D &pt)
{
	CPoint3D::operator=(pt);

	return *this;
}

// get adjacent faces' id of this vertex
bool CVertex::GetAdjacentFacesID(std::vector<int>& adjFaces) const
{
	adjFaces.clear();
	if (m_vecNearEdges.empty())
		return false;

	for (size_t i = 0; i < m_vecNearEdges.size(); ++i)
	{
		CEdge &edge = m_pParent->GetEdgeAt(m_vecNearEdges[i]);

		int fLeft = -1, fRight = -1;
		edge.GetAdjacentFaceId(fLeft, fRight);
		
		//check whether left face is in
		auto iter = std::find(adjFaces.begin(), adjFaces.end(), fLeft);
		if (iter == adjFaces.end())
		{
			adjFaces.push_back(fLeft);
		}

		if (fRight != -1)
		{
			//check whether right face is in
			auto iter = std::find(adjFaces.begin(), adjFaces.end(), fRight);
			if (iter == adjFaces.end())
			{
				adjFaces.push_back(fRight);
			}
		}
	}

	return false;
}

bool CVertex::GetAdjacentVerticesId(std::vector<int> &adjVertices) const
{
	adjVertices.clear();
	adjVertices.resize(m_vecNearEdges.size());
	for (size_t i = 0; i < m_vecNearEdges.size(); ++i)
	{
		CEdge &edge = m_pParent->GetEdgeAt(m_vecNearEdges[i]);
		int id = edge.GetAnotherVertexId(m_id);
		adjVertices.push_back(id);
	}

	return true;
}