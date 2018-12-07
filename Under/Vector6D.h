#pragma once
#include "Matrix.h"
class CVector6D
{
public:
	CVector6D(CVector3f a, CVector3f b){
		m_data1 = a;
		m_data2 = b;
	}
	CVector6D(float a, float b, float c, float d, float e, float f){
		CVector3f temp1(a, b, c);
		CVector3f temp2(d, e, f);
		m_data1 = temp1;
		m_data2 = temp2;
	}

	float operator [](int idx)
	{
		if (idx < 3){
			return m_data1[idx];
		}
		else{
			return m_data2[idx - 3];
		}
	}

	float operator [](int idx) const
	{
		if (idx < 3){
			return m_data1[idx];
		}
		else{
			return m_data2[idx - 3];
		}
	}

	CVector6D operator * (const CVector6D &vec6D) const{
		CVector3f vec3D1 = vec6D.m_data1;
		CVector3f vec3D2 = vec6D.m_data2;
		CVector3f res1 = m_data1%vec3D1;
		CVector3f res2 = m_data2%vec3D2;
		CVector6D res(res1, res2);
		return res;
	}
	CVector3f getData1(){
		return m_data1;
	}
	CVector3f getData2(){
		return m_data2;
	}
private:
	CVector3f m_data1, m_data2;
};
