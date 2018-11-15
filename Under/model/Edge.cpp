#include "stdafx.h"
#include "Edge.h"


CEdge::CEdge(int v0, int v1)
{
	m_vId[0] = v0;
	m_vId[1] = v1;

	m_left = -1;
	m_right = -1;

	m_pParent = nullptr;
}


CEdge::~CEdge()
{
}
