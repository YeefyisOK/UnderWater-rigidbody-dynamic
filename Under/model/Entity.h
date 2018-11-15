#pragma once

#include "Color.h"
#include "BBox.h"
#include "Primitive.h"

#define  POINTCLOUD_TYPE     0x000001
#define  TMESH_TYPE          0x000002
#define  SKELETON_TYPE       0x000004
#define  MESH_TYPE           0x000008
#define  LINESEG_TYPE        0x000010
#define  OBB_TYPE            0x000020
#define  BLENDSHAPE_TYPE     0x000040
#define  ENTITY_TYPE         0xFFFFFF

class CEntity
	:public CPrimitive
{
public:
	CEntity(const std::string &name="");
	~CEntity();

	unsigned int GetObjectType()
	{
		return m_objectType;
	}

	std::string GetName() const
	{
		return m_name;
	}

	void SetName(const std::string &name)
	{
		m_name = name;
	}

	CColor &GetColor()
	{
		return m_color;
	}

	void SetColor(const CColor &color)
	{
		m_color = color;
	}
	void SetColor(byte red, byte green, byte blue)
	{
		m_color.Set(red, green, blue);
	}

	CBBox &GetBBox()
	{
		return m_bbox;
	}

	virtual void CalculateBBox() = 0;
	virtual void Translate(const CVector3D &offset) = 0;

protected:
	std::string m_name;
	CColor m_color;
	CBBox m_bbox;

	unsigned int m_objectType;
};

