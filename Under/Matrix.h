#pragma once
#include "Point3D.h"

class CMatrix
{
public:
	CMatrix(void);
	~CMatrix(void);

	double* operator[](int idx)
	{
		return m_data[idx];
	}

	const double* operator[](int idx) const
	{
		return m_data[idx];
	}

	CMatrix operator % (const CMatrix &matrix) const;//矩阵乘法
	CVector3d operator * (const CVector3d &vector) const;//矩阵乘以向量
	static CMatrix s_GetRotaionMatrix(double angle, const CPoint3D &axis);

private:
	double m_data[3][3];
};

