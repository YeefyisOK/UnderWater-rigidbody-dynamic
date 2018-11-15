#pragma once
class CKirchhoff3D
{
public:
	CKirchhoff3D(double x = 0, double y = 0, double z = 0);
	CKirchhoff3D(double data[3])
	{
		m_data[0] = data[0];
		m_data[1] = data[1];
		m_data[2] = data[2];
	}
	~CKirchhoff3D(void);
private:
	double m_data[3];
};

