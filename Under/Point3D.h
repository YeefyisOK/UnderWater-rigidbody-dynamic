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
	void setData(double data0,double data1,double data2)
	{
		m_data[0] = data0;
		m_data[1] = data1;
		m_data[2] = data2;
	}
	double operator [](int idx)
	{
		return m_data[idx];
	}

	double operator [](int idx) const
	{
		return m_data[idx];
	}

	CPoint3D operator * (const CPoint3D &vec) const;//叉乘

	CPoint3D operator - (const CPoint3D &vec) const;

	CPoint3D& operator += (const CPoint3D &vec);

	double operator % (const CPoint3D &vec) const;//点乘

	CPoint3D& operator /=(const CPoint3D &vec);//向量除以向量 对应相除
	double AngleWith(const CPoint3D &vec) const;

	double Length() const;
	double getm_data(int i)const{
		return m_data[i];
	}
private:
	double m_data[3];
};

typedef CPoint3D CVector3d;

