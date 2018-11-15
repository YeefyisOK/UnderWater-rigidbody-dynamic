#pragma once
#include "Matrix.h"
class CVector6D
{
public:
	CVector6D(CVector3D a, CVector3D b){
		m_data1 = a;
		m_data2 = b;
	}

	double operator [](int idx)
	{
		if (idx < 3){
			return m_data1[idx];
		}
		else{
			return m_data2[idx - 3];
		}
	}

	double operator [](int idx) const
	{
		if (idx < 3){
			return m_data1[idx];
		}
		else{
			return m_data2[idx - 3];
		}
	}

	CVector6D operator * (const CVector6D &vec6D) const{
		CVector3D vec3D1 = vec6D.m_data1;
		CVector3D vec3D2 = vec6D.m_data2;
		CVector3D res1 = m_data1%vec3D1;
		CVector3D res2 = m_data2%vec3D2;
		CVector6D res(res1, res2);
		return res;
	}
	CVector3D getData1(){
		return m_data1;
	}
	CVector3D getData2(){
		return m_data2;
	}
private:
	CVector3D m_data1, m_data2;
};
