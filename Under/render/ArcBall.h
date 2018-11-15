#pragma once
#include "Point3D.h"
#include "Matrix3x3.h"

class CArcBall
{
public:
	CArcBall(int width = 0, int height = 0);
	~CArcBall(void);

	void SetBounds(int width, int height)
	{
		m_nWinWidth = width;
		m_nWinHeight = height;
	}

	CPoint3D WindowToSphere(int x, int y);

	void Click(int x, int y);

	void Drag(int x, int y);

	double * GetRotationData();

private:
	int m_nWinWidth, m_nWinHeight;

	CPoint3D m_ptLast;

	CMatrix3x3 m_matrix;
};

