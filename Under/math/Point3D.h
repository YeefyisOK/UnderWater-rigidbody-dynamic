#pragma once
class CPoint3D
{
public:
	CPoint3D(double x = 0, double y = 0, double z = 0);
	CPoint3D(double data[3])
	{
		m_data[0] = data[0];
		m_data[1] = data[1];
		m_data[2] = data[2];
	}
	~CPoint3D(void);

	double& operator [](int idx)
	{
		return m_data[idx];
	}

	double operator [](int idx) const
	{
		return m_data[idx];
	}

	double X() const
	{
		return m_data[0];
	}

	double Y() const
	{
		return m_data[1];
	}

	double Z() const
	{
		return m_data[1];
	}

	void Set(double x, double y, double z)
	{
		m_data[0] = x;
		m_data[1] = y;
		m_data[2] = z;
	}
	void Set(double coord[3])
	{
		Set(coord[0], coord[1], coord[2]);
	}

	double* Data()
	{
		return m_data;
	}
	const double * Data() const
	{
		return m_data;
	}

	CPoint3D operator * (const CPoint3D &vec) const;
	CPoint3D operator * (double scale) const;

	CPoint3D operator - (const CPoint3D &vec) const;

	CPoint3D operator -() const;

	CPoint3D operator + (const CPoint3D &point) const;

	CPoint3D& operator += (const CPoint3D &vec);
	CPoint3D& operator -= (const CPoint3D &vec);
	CPoint3D& operator *= (double scale);
	CPoint3D& operator /= (double scale)
	{
		return operator *= (1.0 / scale);
	}

	double operator % (const CPoint3D &vec) const;

	double AngleWith(const CPoint3D &vec) const;

	double Length() const;
	double SquareLength() const;

	CPoint3D & Normalize();

	double DistanceWith(const CPoint3D &pt) const
	{
		return sqrt(SquareDistanceWith(pt));
	}

	double SquareDistanceWith(const CPoint3D &pt) const
	{
		return (m_data[0] - pt[0]) * (m_data[0] - pt[0]) +
			(m_data[1] - pt[1]) * (m_data[1] - pt[1]) +
			(m_data[2] - pt[2]) * (m_data[2] - pt[2]);
	}

	void Translate(const CPoint3D &offset)
	{
		m_data[0] += offset[0];
		m_data[1] += offset[1];
		m_data[2] += offset[2];
	}
private:
	double m_data[3];
};

typedef CPoint3D CVector3D;

