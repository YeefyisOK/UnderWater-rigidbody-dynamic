#pragma once
#include "Point3D.h"

class CMatrix
{
public:
	CMatrix(void);
	~CMatrix(void);

	float* operator[](int idx)
	{
		return m_data[idx];
	}

	const float* operator[](int idx) const
	{
		return m_data[idx];
	}

	CMatrix operator % (const CMatrix &matrix) const;//矩阵乘法
	CVector3f operator * (const CVector3f &vector) const;//矩阵乘以向量
	static CMatrix s_GetRotaionMatrix(float angle, const CPoint3D &axis);

private:
	float m_data[3][3];
};

