#pragma once
#include "Matrix.h"
class CVector6D
{
public:
	CVector6D(CVector3d a, CVector3d b){
		m_data1 = a;
		m_data2 = b;
	}
	CVector6D(double a, double b, double c, double d, double e, double f){
		CVector3d temp1(a, b, c);
		CVector3d temp2(d, e, f);
		m_data1 = temp1;
		m_data2 = temp2;
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
		CVector3d vec3D1 = vec6D.m_data1;
		CVector3d vec3D2 = vec6D.m_data2;
		CVector3d res1 = m_data1%vec3D1;
		CVector3d res2 = m_data2%vec3D2;
		CVector6D res(res1, res2);
		return res;
	}
	CVector3d getData1(){
		return m_data1;
	}
	CVector3d getData2(){
		return m_data2;
	}
private:
	CVector3d m_data1, m_data2;
};
