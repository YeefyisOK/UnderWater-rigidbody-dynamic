#pragma once
#include "Point3D.h"

class CObb;

class CMatrix3x3
{
public:

	~CMatrix3x3(void);

	CMatrix3x3()
	{
		memset(m_data[0], 0, 3 * sizeof(double));
		memset(m_data[1], 0, 3 * sizeof(double));
		memset(m_data[2], 0, 3 * sizeof(double));
	}

	CMatrix3x3(const CMatrix3x3 &matrix)
	{
		memcpy(m_data[0], matrix[0], sizeof(double) * 3);
		memcpy(m_data[1], matrix[1], sizeof(double) * 3);
		memcpy(m_data[2], matrix[2], sizeof(double) * 3);
	}


	const CMatrix3x3 & operator = (const CMatrix3x3 &matrix)
	{
		if (this != &matrix)
		{
			memcpy(m_data[0], matrix[0], sizeof(double) * 3);
			memcpy(m_data[1], matrix[1], sizeof(double) * 3);
			memcpy(m_data[2], matrix[2], sizeof(double) * 3);
		}

		return *this;
	}

	CMatrix3x3(double m00, double m01, double m02,
		double m10, double m11, double m12,
		double m20, double m21, double m22)
	{
		m_data[0][0] = m00;   m_data[0][1] = m01;   m_data[0][2] = m02;
		m_data[1][0] = m10;   m_data[1][1] = m11;   m_data[1][2] = m12;
		m_data[2][0] = m20;   m_data[2][1] = m21;   m_data[2][2] = m22;
	}


	double* operator[](int idx)
	{
		return m_data[idx];
	}

	const double* operator[](int idx) const
	{
		return m_data[idx];
	}

	CMatrix3x3 operator * (const CMatrix3x3 &matrix) const;

	CMatrix3x3& operator *= (const CMatrix3x3 &matrix);

	double *RotateData() const;

	void SetIdentity();
	void SetNegativeIdentity();

	void Transpose();

	double Determinant() const;
	CMatrix3x3 Inverse() const;

	bool JacobiEigenv(double eigenValue[3], CMatrix3x3 &eigenVec);

	CMatrix3x3 &FromAxisAngle(CVector3D &axis, double radian)
	{
		*this = s_GetRotaionMatrix(axis, radian);
		return *this;
	}
	CMatrix3x3 &FromAxisAngle(double x, double y, double z, double radian)
	{
		*this = s_GetRotaionMatrix(CVector3D(x, y, z), radian);
		return *this;
	}

	void FromTwoAxis(const CVector3D &from, const CVector3D &to);

	

	//conversion to euler angle
	bool ToEulerAnglesXYZ(double& xAngle, double& yAngle, double& zAngle) const;

	bool ToEulerAnglesXZY(double& xAngle, double& zAngle, double& yAngle) const;

	bool ToEulerAnglesYXZ(double& yAngle, double& xAngle, double& zAngle) const;

	bool ToEulerAnglesYZX(double& yAngle, double& zAngle, double& xAngle) const;

	bool ToEulerAnglesZXY(double& zAngle, double& xAngle, double& yAngle) const;

	bool ToEulerAnglesZYX(double& zAngle, double& yAngle, double& xAngle) const;

	//conversion from euler angle
	void FromEulerAnglesXYZ(const double& fYAngle, const double& fPAngle, const double& fRAngle);

	void FromEulerAnglesXZY(const double& fYAngle, const double& fPAngle, const double& fRAngle);

	void FromEulerAnglesYXZ(const double& fYAngle, const double& fPAngle, const double& fRAngle);

	void FromEulerAnglesYZX(const double& fYAngle, const double& fPAngle, const double& fRAngle);

	void FromEulerAnglesZXY(const double& fYAngle, const double& fPAngle, const double& fRAngle);

	void FromEulerAnglesZYX(const double& fYAngle, const double& fPAngle, const double& fRAngle);


	friend CVector3D operator * (const CMatrix3x3 &mtx, const CVector3D &vec);

	static CMatrix3x3 s_GetRotaionMatrix( const CPoint3D &axis, double radian);

	static CMatrix3x3 s_GetCovMatrix(const std::vector<CPoint3D> &vecPnts);

private:
	double m_data[3][3];
};

