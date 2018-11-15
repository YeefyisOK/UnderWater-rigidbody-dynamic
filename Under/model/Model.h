#pragma once

#include "BBox.h"
#include "Entity.h"
#include "Mesh.h"


class CModel
{
public:
	CModel();
	~CModel();

	static CModel& GetInstance()
	{
		return m_model;
	}

	void clear()
	{
		for (size_t i = 0; i < m_vecEnities.size(); ++i)
		{
			delete m_vecEnities[i];
		}
		m_vecEnities.clear();

		m_bbox.SetEmpty();
	}

	void AddEntity(CEntity *pEntity, bool bCalculateBBox = true);
	

	void RemoveEntity(CEntity *pEntity);
	
	size_t GetNumOfEntities() const 
	{
		return m_vecEnities.size();
	}

	CEntity *GetEntityAt(int idx)
	{
		return m_vecEnities[idx];
	}

	

	const CBBox &GetBBox() const
	{
		return m_bbox;
	}

	void CalculateBBox(bool bCalcEntity = true);

	
	CMesh *GetMesh();
	
	void ClearSelectedTag();
private:
	std::vector<CEntity *> m_vecEnities;

	CBBox m_bbox;

	static CModel m_model;
};

