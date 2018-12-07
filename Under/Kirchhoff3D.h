#pragma once
class CKirchhoff3D
{
public:
	CKirchhoff3D(float x = 0, float y = 0, float z = 0);
	CKirchhoff3D(float data[3])
	{
		m_data[0] = data[0];
		m_data[1] = data[1];
		m_data[2] = data[2];
	}
	~CKirchhoff3D(void);
private:
	float m_data[3];
};

