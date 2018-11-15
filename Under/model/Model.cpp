#include "stdafx.h"
#include "Model.h"
#include "Entity.h"

CModel CModel::m_model;

CModel::CModel()
{
}


CModel::~CModel()
{
	for (size_t i = 0; i < m_vecEnities.size(); ++i)
	{
		delete m_vecEnities[i];
	}
	m_vecEnities.clear();
}

void CModel::AddEntity(CEntity *pEntity, bool bCalculateBBox /* = true */)
{
	if (pEntity)
	{
		if (bCalculateBBox)
		{
			pEntity->CalculateBBox();
		}
		
		if (m_vecEnities.empty())
		{
			m_bbox = pEntity->GetBBox();
		} 
		else
		{
			m_bbox.MergeWith(pEntity->GetBBox());
		}
		m_vecEnities.push_back(pEntity);
	}
}

void CModel::RemoveEntity(CEntity *pEntity)
{
	if (pEntity)
	{
		auto iter = std::find(m_vecEnities.begin(), m_vecEnities.end(), pEntity);
		m_vecEnities.erase(iter);
		delete pEntity;
	}
}


CMesh *CModel::GetMesh()
{
	for (size_t i = 0; i < m_vecEnities.size(); ++i)
	{
		if (m_vecEnities[i]->GetObjectType() == MESH_TYPE)
		{
			return (CMesh *)m_vecEnities[i];
		}
	}

	return nullptr;
}


void CModel::CalculateBBox(bool bCalcEntity /* = true */)
{
	if (m_vecEnities.empty())
	{
		return;
	}

	if (bCalcEntity)
	{
		m_vecEnities[0]->CalculateBBox();
	}
	m_bbox = m_vecEnities[0]->GetBBox();
	for (size_t i = 1; i < m_vecEnities.size(); ++i)
	{
		if (bCalcEntity)
		{
			m_vecEnities[i]->CalculateBBox();
		}
		m_bbox.MergeWith(m_vecEnities[i]->GetBBox());
	}
}

void CModel::ClearSelectedTag()
{
	for (size_t i = 0; i < m_vecEnities.size(); ++i)
	{
		m_vecEnities[i]->SetSelected(false);
	}
}