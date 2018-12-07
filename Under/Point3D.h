#pragma once
class CPoint3D
{
public:
	CPoint3D(float x = 0, float y = 0, float z = 0);
	CPoint3D(float data[3])
	{
		m_data[0] = data[0];
		m_data[1] = data[1];
		m_data[2] = data[2];
	}
	~CPoint3D(void);
	void setData(float data0,float data1,float data2)
	{
		m_data[0] = data0;
		m_data[1] = data1;
		m_data[2] = data2;
	}
	float operator [](int idx)
	{
		return m_data[idx];
	}

	float operator [](int idx) const
	{
		return m_data[idx];
	}

	CPoint3D operator * (const CPoint3D &vec) const;//叉乘

	CPoint3D operator - (const CPoint3D &vec) const;

	CPoint3D& operator += (const CPoint3D &vec);

	float operator % (const CPoint3D &vec) const;//点乘

	CPoint3D& operator /=(const CPoint3D &vec);//向量除以向量 对应相除
	float AngleWith(const CPoint3D &vec) const;

	float Length() const;
	float getm_data(int i)const{
		return m_data[i];
	}
private:
	float m_data[3];
};

typedef CPoint3D CVector3f;

