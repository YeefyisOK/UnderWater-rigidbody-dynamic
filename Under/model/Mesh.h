#pragma once
#include "Entity.h"

#include "Vertex.h"
#include "Edge.h"
#include "Face.h"

class CMesh
	:public CEntity
{
public:
	CMesh();
	~CMesh();

	CMesh(const CMesh &mesh);
	CMesh & operator = (const CMesh &mesh);

	size_t GetNumOfVertices()
	{
		return m_vertices.size();
	}
	void SetNumOfVertices(size_t num)
	{
		m_vertices.resize(num);
	}

	size_t GetNumOfEdges()
	{
		return m_vecEdges.size();
	}
	void SetNumOfEdges(size_t num)
	{
		m_vecEdges.resize(num);
	}

	size_t GetNumOfFaces()
	{
		return m_vecFaces.size();
	}
	void SetNumOfFaces(size_t num)
	{
		m_vecFaces.resize(num);
	}

	CVertex& GetVertexAt(size_t idx)
	{
		return m_vertices[idx];
	}
	void SetVertexAt(int idx, double x, double y, double z)
	{
		m_vertices[idx].SetId(idx);
		m_vertices[idx].SetParent(this);
		m_vertices[idx].Set(x, y, z);
	}
	void SetVertexAt(int idx, double coord[3])
	{
		SetVertexAt(idx, coord[0], coord[1], coord[2]);
	}

	void SetVertexAt(int idx, const CPoint3D &pt)
	{
		SetVertexAt(idx, pt[0], pt[1], pt[2]);
	}

	std::vector<CVertex> &GetAllVertices()
	{
		return m_vertices;
	}

	CEdge& GetEdgeAt(size_t idx)
	{
		return m_vecEdges[idx];
	}
	int AddEdge(int vId0, int vId1)
	{
		int id = (int)m_vecEdges.size();

		m_vecEdges.push_back(CEdge(vId0, vId1));
		
		m_vecEdges.back().SetId(id);
		m_vecEdges.back().SetParent(this);

		return id;
	}
	std::vector<CEdge>& GetAllEdges()
	{
		return m_vecEdges;
	}

	CFace &GetFaceAt(size_t idx)
	{
		return m_vecFaces[idx];
	}
	void SetFaceAt(int idx, int v0, int v1, int v2)
	{
		m_vecFaces[idx].SetId(idx);
		m_vecFaces[idx].SetParent(this);
		m_vecFaces[idx].SetVerticeId(v0, v1, v2);
	}
	void SetFaceAt(int idx, int vId[3])
	{
		SetFaceAt(idx, vId[0], vId[1], vId[2]);
	}
	std::vector<CFace> &GetAllFaces()
	{
		return m_vecFaces;
	}

	void SetTexture(cv::Mat image)
	{
		m_image = image;
	}

	void CalculateBBox();

	void Translate(const CVector3D &offset);

	void Scale(const CPoint3D &ptCenter, double scale);

	bool HasTexture() const
	{
		return !m_image.empty();
	}
	bool IsTextureIntited() const
	{
		return m_textureId > 0;
	}
	unsigned int GetTextureId() const
	{
		return m_textureId;
	}

	void SetTextureId(unsigned int id)
	{
		m_textureId = id;
	}

	cv::Mat &GetTexture()
	{
		return m_image;
	}

	int GetEdgeIDFromVertices(int v0, int v1);

	void ComputeFaceNormal();
	void ComputeFaceVertexNormal();

	void CalculateEdgesFromFacesAndVertices();

private:
	std::vector<CVertex> m_vertices;
	std::vector<CEdge> m_vecEdges;
	std::vector<CFace > m_vecFaces;

	cv::Mat m_image;
	unsigned int m_textureId;
};

