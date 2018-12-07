#include "Point3D.h"
#include <math.h>

CPoint3D::CPoint3D(float x, float y, float z)
{
	m_data[0] = x;
	m_data[1] = y;
	m_data[2] = z;
}


CPoint3D::~CPoint3D(void)
{
}


CPoint3D CPoint3D::operator *(const CPoint3D &vec) const
{
	float x = m_data[1] * vec.m_data[2] - vec.m_data[1] * m_data[2];
	float y = vec.m_data[0] * m_data[2] - m_data[0] * vec.m_data[2];
	float z = m_data[0] * vec.m_data[1] - vec.m_data[0] * m_data[1];

	return CPoint3D(x, y, z);
}

CPoint3D CPoint3D::operator -(const CPoint3D &point) const
{
	return CPoint3D(m_data[0] - point[0], m_data[1] - point[1], m_data[2] - point[2]);
}

CPoint3D& CPoint3D::operator +=(const CPoint3D &point)
{
	m_data[0] += point[0];
	m_data[1] += point[1];
	m_data[2] += point[2];

	return *this;
}

float CPoint3D::operator %(const CPoint3D &vec) const
{
	float dot = m_data[0] * vec.m_data[0] + m_data[1] * vec.m_data[1] + m_data[2] * vec.m_data[2];
	return dot;
}

CPoint3D& CPoint3D::operator /= (const CPoint3D &vec){//向量除以向量 对应相除
	m_data[0] /= vec.getm_data(0);
	m_data[1] /= vec.getm_data(1);
	m_data[2] /= vec.getm_data(2);

	return *this;
}
float CPoint3D::AngleWith(const CPoint3D &vec) const
{
	CPoint3D cross = (*this) * vec;
	float len = cross.Length();

	float dot = (*this) % vec;

	float angle = atan2(len, dot);

	return angle;
}

float CPoint3D::Length() const
{
	float len = m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2];
	return sqrt(len);

}