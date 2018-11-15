#include "StdAfx.h"
#include "ArcBall.h"
#include <math.h>

CArcBall::CArcBall(int width, int height)
{
	m_nWinWidth = width;
	m_nWinHeight = height;

	m_matrix.SetIdentity();
}


CArcBall::~CArcBall(void)
{
}

CPoint3D CArcBall::WindowToSphere(int x, int y)
{
	double xScale = 2.0 / m_nWinWidth;
	double yScale = 2.0 / m_nWinHeight;

	double xPos = x * xScale - 1;
	double yPos = 1 - y * yScale;
	double zPos = 0.0;

	double len = xPos * xPos + yPos * yPos;
	if (len < 1.0)
	{
		zPos = sqrt(1.0 - len);
	}

	return CPoint3D(xPos, yPos, zPos);
}


void CArcBall::Click(int x, int y)
{
	m_ptLast = WindowToSphere(x, y);
}

void CArcBall::Drag(int x, int y)
{
	CPoint3D ptCur = WindowToSphere(x, y);

	CPoint3D axis = m_ptLast * ptCur;
	if (axis.Length() < 1.0e-7)
	{
		return;
	}

	double angle = 2.0 * m_ptLast.AngleWith(ptCur);

	//matrix
	CMatrix3x3 matrix = CMatrix3x3::s_GetRotaionMatrix(axis, angle);

	m_matrix = matrix * m_matrix;

	m_ptLast = ptCur;
}

double *CArcBall::GetRotationData()
{
	static double data[16];

	data[0] = m_matrix[0][0];
	data[1] = m_matrix[1][0];
	data[2] = m_matrix[2][0];
	data[3] = 0;

	data[4] = m_matrix[0][1];
	data[5] = m_matrix[1][1];
	data[6] = m_matrix[2][1];
	data[7] = 0;

	data[8] = m_matrix[0][2];
	data[9] = m_matrix[1][2];
	data[10] = m_matrix[2][2];
	data[11] = 0;

	data[12] = 0;
	data[13] = 0;
	data[14] = 0;
	data[15] = 1;

	return data;
}