#include "stdafx.h"
#include "Mesh.h"


CMesh::CMesh()
{
	m_objectType = MESH_TYPE;
	m_textureId = 0;

	//m_color.Set(128, 0, 255);
	m_color.Set(192, 192, 192);
}


CMesh::~CMesh()
{
}

CMesh::CMesh(const CMesh &mesh)
	:CEntity(mesh)
{
	m_vertices = mesh.m_vertices;
	m_vecEdges = mesh.m_vecEdges;
	m_vecFaces = mesh.m_vecFaces;

	for (size_t i = 0; i < m_vecFaces.size(); ++i)
	{
		m_vecFaces[i].SetParent(this);
	}

	for (size_t i = 0; i < m_vecEdges.size(); ++i)
	{
		m_vecEdges[i].SetParent(this);
	}

	for (size_t i = 0; i < m_vertices.size(); ++i)
	{
		m_vertices[i].SetParent(this);
	}

	m_objectType = MESH_TYPE;
}

CMesh &CMesh::operator=(const CMesh &mesh)
{
	if (this != &mesh)
	{
		CEntity::operator = (mesh);

		m_vertices = mesh.m_vertices;
		m_vecEdges = mesh.m_vecEdges;
		m_vecFaces = mesh.m_vecFaces;

		for (size_t i = 0; i < m_vecFaces.size(); ++i)
		{
			m_vecFaces[i].SetParent(this);
		}

		for (size_t i = 0; i < m_vecEdges.size(); ++i)
		{
			m_vecEdges[i].SetParent(this);
		}

		for (size_t i = 0; i < m_vertices.size(); ++i)
		{
			m_vertices[i].SetParent(this);
		}

		m_objectType = MESH_TYPE;
	}

	return *this;
}

void CMesh::CalculateBBox()
{
	if (m_vertices.empty())
	{
		return;
	}

	CPoint3D ptMin = m_vertices[0];
	CPoint3D ptMax = ptMin;

	for (size_t i = 1; i < m_vertices.size(); ++i)
	{
		//x
		if (m_vertices[i][0] < ptMin[0])
		{
			ptMin[0] = m_vertices[i][0];
		}
		else if (m_vertices[i][0] > ptMax[0])
		{
			ptMax[0] = m_vertices[i][0];
		}

		//y
		if (m_vertices[i][1] < ptMin[1])
		{
			ptMin[1] = m_vertices[i][1];
		}
		else if (m_vertices[i][1] > ptMax[1])
		{
			ptMax[1] = m_vertices[i][1];
		}

		//z
		if (m_vertices[i][2] < ptMin[2])
		{
			ptMin[2] = m_vertices[i][2];
		}
		else if (m_vertices[i][2] > ptMax[2])
		{
			ptMax[2] = m_vertices[i][2];
		}
	}

	m_bbox.Set(ptMin, ptMax);
}

void CMesh::Translate(const CVector3D &offset)
{
	for (size_t i = 0; i < m_vertices.size(); ++i)
	{
		m_vertices[i].Translate(offset);
	}
	m_bbox.Translate(offset);
}

void CMesh::ComputeFaceNormal()
{
	for (size_t i = 0; i < m_vecFaces.size(); ++i)
	{
		m_vecFaces[i].ComputeNormal();
	}
}

void CMesh::ComputeFaceVertexNormal()
{
	ComputeFaceNormal();

	std::vector<int> adjFacesId;
	for (size_t i = 0; i < m_vertices.size(); ++i)
	{
		//1-ring
		m_vertices[i].GetAdjacentFacesID(adjFacesId);

		CVector3D normal;
		for (size_t j = 0; j < adjFacesId.size(); ++j)
		{
			normal += m_vecFaces[adjFacesId[j]].GetNormal();
		}
		m_vertices[i].SetNormal(normal);
	}
}

void CMesh::CalculateEdgesFromFacesAndVertices()
{
	int vNext = -1;
	int eId[3];
	for (size_t i = 0; i < m_vecFaces.size(); ++i)
	{
		int *vId = m_vecFaces[i].GetVerticeId();

		for (int j = 0; j < 3; ++j)
		{
			vNext = vId[(j + 1) % 3];
			eId[j] = GetEdgeIDFromVertices(vId[j], vNext);
			if (eId[j] == -1)
			{
				eId[j] = AddEdge(vId[j], vNext);
				m_vertices[vId[j]].AddNearEdge(eId[j]);
				m_vertices[vNext].AddNearEdge(eId[j]);
			}
		}

		m_vecFaces[i].SetEdgesId(eId);
	}

	for (size_t i = 0; i < m_vecFaces.size(); ++i)
	{
		int *eId = m_vecFaces[i].GetEdgesId();
		for (int j = 0; j < 3; ++j)
		{
			CEdge &edge = m_vecEdges[eId[j]];
			edge.AddFace(i);
			
		}
	}
}

int CMesh::GetEdgeIDFromVertices(int v0, int v1)
{
	std::vector<int> & vecNearEdges = m_vertices[v0].GetNearEdges();
	for (size_t i = 0; i < vecNearEdges.size(); i++)
	{
		CEdge &edge = m_vecEdges[vecNearEdges[i]];
		int v = edge.GetAnotherVertexId(v0);

		if (v == v1)
		{
			return vecNearEdges[i];
		}
	}

	return -1;
}

void CMesh::Scale(const CPoint3D &ptCenter, double scale)
{
	for (size_t i = 0; i < m_vertices.size(); ++i)
	{
		m_vertices[i] = (m_vertices[i] - ptCenter) * scale + ptCenter;
	}

	//m_bbox.Scale(ptCenter, scale);
	CalculateBBox();
}