#include "StdAfx.h"
#include "Point3D.h"
#include <math.h>

CPoint3D::CPoint3D(double x, double y, double z)
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
	double x = m_data[1] * vec.m_data[2] - vec.m_data[1] * m_data[2];
	double y = vec.m_data[0] * m_data[2] - m_data[0] * vec.m_data[2];
	double z = m_data[0] * vec.m_data[1] - vec.m_data[0] * m_data[1];

	return CPoint3D(x, y, z);
}

CPoint3D CPoint3D::operator* (double scale) const
{
	return CPoint3D(m_data[0] * scale, m_data[1] * scale, m_data[2] * scale);
}

CPoint3D CPoint3D::operator -(const CPoint3D &point) const
{
	return CPoint3D(m_data[0] - point[0], m_data[1] - point[1], m_data[2] - point[2]);
}

CPoint3D CPoint3D::operator-() const
{
	return CPoint3D(-m_data[0], -m_data[1], -m_data[2]);
}

CPoint3D CPoint3D::operator +(const CPoint3D &point) const
{
	return CPoint3D(m_data[0] + point[0], m_data[1] + point[1], m_data[2] + point[2]);
}

CPoint3D& CPoint3D::operator +=(const CPoint3D &point)
{
	m_data[0] += point[0];
	m_data[1] += point[1];
	m_data[2] += point[2];

	return *this;
}

CPoint3D& CPoint3D::operator -=(const CPoint3D &point)
{
	m_data[0] -= point[0];
	m_data[1] -= point[1];
	m_data[2] -= point[2];

	return *this;
}

CPoint3D & CPoint3D::operator*= (double scale)
{
	m_data[0] *= scale;
	m_data[1] *= scale;
	m_data[2] *= scale;

	return *this;
}

double CPoint3D::operator %(const CPoint3D &vec) const
{
	double dot = m_data[0] * vec.m_data[0] + m_data[1] * vec.m_data[1] + m_data[2] * vec.m_data[2];
	return dot;
}

double CPoint3D::AngleWith(const CPoint3D &vec) const
{
	CPoint3D cross = (*this) * vec;
	double len = cross.Length();

	double dot = (*this) % vec;

	double angle = atan2(len, dot);

	return angle;
}

double CPoint3D::SquareLength() const
{
	double len = m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2];
	return len;
}

double CPoint3D::Length() const
{
	double len = m_data[0] * m_data[0] + m_data[1] * m_data[1] + m_data[2] * m_data[2];
	return sqrt(len);

}

CPoint3D & CPoint3D::Normalize()
{
	double len = Length();
	if (len < 1.0e-8)
	{
		return *this;
	}

	len = 1.0 / len;
	m_data[0] *= len;
	m_data[1] *= len;
	m_data[2] *= len;

	return *this;
}