#pragma once
#include "Primitive.h"

class CMesh;

class CTriangle
	:public CPrimitive
{
public:
	CTriangle();
	~CTriangle();

private:
	int m_vId[3];

	CMesh *m_pParent;
};

