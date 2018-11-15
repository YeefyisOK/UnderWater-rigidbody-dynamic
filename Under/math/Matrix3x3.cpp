#include "StdAfx.h"
#include "Matrix3x3.h"
#include <math.h>


CMatrix3x3::~CMatrix3x3(void)
{
}

void CMatrix3x3::SetIdentity()
{
	memset(m_data[0], 0, sizeof(double) * 9);
	m_data[0][0] = 1.0;
	m_data[1][1] = 1.0;
	m_data[2][2] = 1.0;
}

void CMatrix3x3::SetNegativeIdentity()
{
	memset(m_data[0], 0, sizeof(double) * 9);
	m_data[0][0] = -1.0;
	m_data[1][1] = -1.0;
	m_data[2][2] = -1.0;
}

CMatrix3x3 CMatrix3x3::operator *(const CMatrix3x3 &matrix) const
{
#define MATRIX_MULT(i, j)  m_data[i][0] * matrix[0][j] + m_data[i][1] * matrix[1][j] + m_data[i][2] * matrix[2][j]
	CMatrix3x3 result;
	result[0][0] = MATRIX_MULT(0, 0);
	result[0][1] = MATRIX_MULT(0, 1);
	result[0][2] = MATRIX_MULT(0, 2);

	result[1][0] = MATRIX_MULT(1, 0);
	result[1][1] = MATRIX_MULT(1, 1);
	result[1][2] = MATRIX_MULT(1, 2);

	result[2][0] = MATRIX_MULT(2, 0);
	result[2][1] = MATRIX_MULT(2, 1);
	result[2][2] = MATRIX_MULT(2, 2);


	return result;
}

CMatrix3x3& CMatrix3x3::operator*=(const CMatrix3x3 &matrix)
{
#define MATRIX_MULT(i, j)  data[0] * matrix[0][j] +data[1] * matrix[1][j] + data[2] * matrix[2][j]
	double data[3];

	data[0] = m_data[0][0]; data[1] = m_data[0][1]; data[2] = m_data[0][2];
	m_data[0][0] = MATRIX_MULT(0, 0);
	m_data[0][1] = MATRIX_MULT(0, 1);
	m_data[0][2] = MATRIX_MULT(0, 2);

	data[0] = m_data[1][0]; data[1] = m_data[1][1]; data[2] = m_data[1][2];
	m_data[1][0] = MATRIX_MULT(1, 0);
	m_data[1][1] = MATRIX_MULT(1, 1);
	m_data[1][2] = MATRIX_MULT(1, 2);

	data[0] = m_data[2][0]; data[1] = m_data[2][1]; data[2] = m_data[2][2];
	m_data[2][0] = MATRIX_MULT(2, 0);
	m_data[2][1] = MATRIX_MULT(2, 1);
	m_data[2][2] = MATRIX_MULT(2, 2);

	return *this;
}

void CMatrix3x3::Transpose()
{
	std::swap(m_data[0][1], m_data[1][0]);
	std::swap(m_data[0][2], m_data[2][0]);
	std::swap(m_data[1][2], m_data[2][1]);
}

double *CMatrix3x3::RotateData() const
{
	static double data[16];

	data[0] = m_data[0][0];
	data[1] = m_data[1][0];
	data[2] = m_data[2][0];
	data[3] = 0;

	data[4] = m_data[0][1];
	data[5] = m_data[1][1];
	data[6] = m_data[2][1];
	data[7] = 0;

	data[8] = m_data[0][2];
	data[9] = m_data[1][2];
	data[10] = m_data[2][2];
	data[11] = 0;

	data[12] = 0;
	data[13] = 0;
	data[14] = 0;
	data[15] = 1;

	return data;
}

// bool CMatrix3x3::JacobiEigenv(double eigenValue[3], CMatrix3x3 &eigenVec)
// {
// 	Eigen::Matrix3d temp;
// 	temp(0, 0) = m_data[0][0];
// 	temp(0, 1) = m_data[0][1];
// 	temp(0, 2) = m_data[0][2];
// 	temp(1, 0) = m_data[1][0];
// 	temp(1, 1) = m_data[1][1];
// 	temp(1, 2) = m_data[1][2];
// 	temp(2, 0) = m_data[2][0];
// 	temp(2, 1) = m_data[2][1];
// 	temp(2, 2) = m_data[2][2];
// 
// //	temp << 1.23, 2.12, -4.2, 2.12, -5.6, 8.79, -4.2, 8.79, 7.3;
// 
// 
// 	Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigenSolver(temp);
// 	if (eigenSolver.info() == Eigen::Success) {
// 
// 		//sort
// 		auto &value = eigenSolver.eigenvalues();
// 
// 		std::vector<std::pair<double, int> > vecValuePairs;
// 		vecValuePairs.push_back(std::make_pair(fabs(value(0, 0)), 0));
// 		vecValuePairs.push_back(std::make_pair(fabs(value(1, 0)), 1));
// 		vecValuePairs.push_back(std::make_pair(fabs(value(2, 0)), 2));
// 
// 		std::sort(vecValuePairs.begin(), vecValuePairs.end());
// 
// 		eigenValue[0] = value(vecValuePairs[2].second, 0);
// 		eigenValue[1] = value(vecValuePairs[1].second, 0);
// 		eigenValue[2] = value(vecValuePairs[0].second, 0);
// 
// 		auto &vectors = eigenSolver.eigenvectors();
// 		
// 		eigenVec[0][0] = vectors(0, vecValuePairs[2].second);
// 		eigenVec[1][0] = vectors(1, vecValuePairs[2].second);
// 		eigenVec[2][0] = vectors(2, vecValuePairs[2].second);
// 
// 		eigenVec[0][1] = vectors(0, vecValuePairs[1].second);
// 		eigenVec[1][1] = vectors(1, vecValuePairs[1].second);
// 		eigenVec[2][1] = vectors(2, vecValuePairs[1].second);
// 
// 		eigenVec[0][2] = vectors(0, vecValuePairs[0].second);
// 		eigenVec[1][2] = vectors(1, vecValuePairs[0].second);
// 		eigenVec[2][2] = vectors(2, vecValuePairs[0].second);
// 
// 		return true;
// 	}
// 
// 	return false;
// }

void CMatrix3x3::FromTwoAxis(const CVector3D &from, const CVector3D &to)
{
	CVector3D fromAxis = from;
	CVector3D toAxis = to;

	fromAxis.Normalize();
	toAxis.Normalize();

	// dot product
	double cosAngle = fromAxis % toAxis;
	//if angle is too little then return 3X3 identity matrix
	if (cosAngle > 1.0 - FLT_EPSILON)
	{
		SetIdentity();
		return;
	}
	else if (cosAngle < -1.0 + FLT_EPSILON)
	{
		//negative identity
		SetNegativeIdentity();    //not a rotation matrix, how to improve it ?
		return;
	}

	// cross product
	CVector3D utVec = fromAxis * toAxis;

	//normalize the rotation axis	
	// 	double mag2 = utVec.SquareLength();
	// 	if (mag2 < CConfig::m_zeroTolerance)
	// 	{
	// 		SetIdentity();
	// 		return;
	// 	}
	// 
	// 	double s = sqrt(mag2);
	// 
	// 	double mag = 1 / s;
	// 
	// 	double x = utVec.X() * mag;
	// 	double y = utVec.Y() * mag;
	// 	double z = utVec.Z() * mag;

	utVec.Normalize();
	double x = utVec.X();
	double y = utVec.Y();
	double z = utVec.Z();

	double c = cosAngle;
	double s = sqrt(1.0 - c * c);
	double t = 1.0 - c;

	double tx = t * x;
	double ty = t * y;
	double tz = t * z;
	double sx = s * x;
	double sy = s * y;
	double sz = s * z;

	//-----------------------------------------------------------
	//		| t*x*x + c		t*x*y - s*z		t*x*z + s*y |
	//		|											|
	//	R = | t*x*y + s*z	t*y*y + c		t*y*z - s*x |
	//		|											|
	//		| t*x*z - s*y	t*y*z + s*x		t*z*z + c	|
	//
	// where c = cos(theta), s = sin(theta), t = 1 - c and(x, y, z) is a unit
	// vector on the axis of rotation.
	//-----------------------------------------------------------

	// row one
	m_data[0][0] = tx * x + c;
	m_data[0][1] = tx * y - sz;
	m_data[0][2] = tx * z + sy;

	// row two
	m_data[1][0] = m_data[0][1] + sz + sz;	// tx * y + sz
	m_data[1][1] = ty * y + c;
	m_data[1][2] = ty * z - sx;

	// row three
	m_data[2][0] = m_data[0][2] - sy - sy;	// tx * z + sy
	m_data[2][1] = m_data[1][2] + sx + sx;	// ty * z + sx
	m_data[2][2] = tz * z + c;

	//double det = Determinant();
	//assert(CConfig::IsZero(Determinant() - 1));
}

CMatrix3x3 CMatrix3x3::s_GetRotaionMatrix( const CPoint3D &axis, double radian)
{
	double x = axis[0];
	double y = axis[1];
	double z = axis[2];

#define  ZERO_TOL    1.0e-7
	double s = sin(radian);
	double c = cos(radian);

	CMatrix3x3 matrix;
	if (fabs(z) < ZERO_TOL)
	{
		if (fabs(y) < ZERO_TOL)
		{
			//if the axis is (0, 0, 0), assume to be identity matrix
			if(fabs(x) < ZERO_TOL)   
				return matrix;

			//rotation axis is (1, 0, 0)
			matrix[0][0] = 1.0;
			matrix[1][1] = c;  matrix[2][2] = c;

			if(x > 0)
			{
				matrix[1][2] = -s; matrix[2][1] = s;
			}
			else
			{
				matrix[1][2] = s; matrix[2][1] = -s;
			}

			return matrix;

		}
		else if (fabs(x) < ZERO_TOL)
		{
			//rotation axis is (0, 1, 0)
			matrix[1][1] = 1.0;
			matrix[0][0] = c;  matrix[2][2] = c;

			if(y > 0)
			{
				matrix[0][2] = s;  matrix[2][0] = -s;
			}
			else
			{
				matrix[0][2] = -s;  matrix[2][0] = s;
			}			

			return matrix;
		}
	}
	else if (fabs(y) < ZERO_TOL)
	{
		if (fabs(x) < ZERO_TOL)
		{
			//rotation axis is (0, 0, 1)
			matrix[2][2] = 1.0;
			matrix[0][0] = c;  matrix[1][1] = c;

			if (z > 0)
			{
				matrix[0][1] = -s; matrix[1][0] = s;
			}
			else
			{
				matrix[0][1] = s; matrix[1][0] = -s;
			}			

			return matrix;
		}
	}

	//common case

	//normalize the rotation axis
	double mag = sqrt(x * x + y * y + z * z);
	mag = 1.0 / mag;
	x *= mag;
	y *= mag;
	z *= mag;

	double t = 1.0 - c;

	double tx = t * x;
	double ty = t * y;
	double tz = t * z;
	double sx = s * x;
	double sy = s * y;
	double sz = s * z;

	//-----------------------------------------------------------
	//		| t*x*x + c		t*x*y - s*z		t*x*z + s*y |
	//		|											|
	//	R = | t*x*y + s*z	t*y*y + c		t*y*z - s*x |
	//		|											|
	//		| t*x*z - s*y	t*y*z + s*x		t*z*z + c	|
	//
	// where c = cos(theta), s = sin(theta), t = 1 - c and(x, y, z) is a unit
	// vector on the axis of rotation.
	//-----------------------------------------------------------

	// row one
	matrix[0][0] = tx * x + c;
	matrix[0][1] = tx * y - sz;
	matrix[0][2] = tx * z + sy;

	// row two
	matrix[1][0] = matrix[0][1] + sz + sz;	// tx * y + sz
	matrix[1][1] = ty * y + c;
	matrix[1][2] = ty * z - sx;

	// row three
	matrix[2][0] = matrix[0][2] - sy - sy;	// tx * z - sy
	matrix[2][1] = matrix[1][2] + sx + sx;	// ty * z + sx
	matrix[2][2] = tz * z + c;

	return matrix;

}

//-----------------------------------------------------------------------
bool CMatrix3x3::ToEulerAnglesXYZ(double& xAngle, double& yAngle,
	double& zAngle) const
{
	// rot =  cy*cz          -cy*sz           sy
	//        cz*sx*sy+cx*sz  cx*cz-sx*sy*sz -cy*sx
	//       -cx*cz*sy+sx*sz  cz*sx+cx*sy*sz  cx*cy

	if (m_data[0][2] > 1.0 - FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		yAngle = HALF_PI;
		zAngle = 0.0;  // any angle works
		double fRpY = atan2(m_data[1][0], m_data[1][1]);
		xAngle = fRpY - zAngle;
		return false;
	}
	else if (m_data[0][2] < -1.0 + FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		yAngle = -HALF_PI;
		zAngle = 0.0;      // any angle works
		double fRmY = atan2(m_data[1][0], m_data[1][1]);
		xAngle = zAngle - fRmY;
		return false;
	}

	yAngle = double(asin(m_data[0][2]));
	xAngle = atan2(-m_data[1][2], m_data[2][2]);
	zAngle = atan2(-m_data[0][1], m_data[0][0]);

	return true;
}
//-----------------------------------------------------------------------
bool CMatrix3x3::ToEulerAnglesXZY(double& xAngle, double& zAngle,
	double& yAngle) const
{
	// rot =  cy*cz          -sz              cz*sy
	//        sx*sy+cx*cy*sz  cx*cz          -cy*sx+cx*sy*sz
	//       -cx*sy+cy*sx*sz  cz*sx           cx*cy+sx*sy*sz

	double sz = -m_data[0][1];

	if (sz > 1.0 - FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		zAngle = HALF_PI;
		yAngle = 0.0;  // any angle works

		double fRpY = atan2(-m_data[1][0], m_data[1][2]);
		xAngle = fRpY - yAngle;
		return false;
	}
	else if (sz < -1.0 + FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		zAngle = -HALF_PI;
		yAngle = 0.0;  // any angle works

		double fRmY = atan2(-m_data[1][0], m_data[1][2]);
		xAngle = yAngle - fRmY;
		return false;
	}

	zAngle = asin(sz);
	xAngle = atan2(m_data[2][1], m_data[1][1]);
	yAngle = atan2(m_data[0][2], m_data[0][0]);
	return true;
}

//-----------------------------------------------------------------------
bool CMatrix3x3::ToEulerAnglesYXZ(double& yAngle, double& xAngle,
	double& zAngle) const
{
	// rot =  cy*cz+sx*sy*sz  cz*sx*sy-cy*sz  cx*sy
	//        cx*sz           cx*cz          -sx
	//       -cz*sy+cy*sx*sz  cy*cz*sx+sy*sz  cx*cy

	double sx = -m_data[1][2];
	if (sx > 1.0 - FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		xAngle = HALF_PI;
		zAngle = 0.0;  // any angle works

		double fRpY = atan2(-m_data[0][1], m_data[0][0]);
		yAngle = fRpY - zAngle;
		return false;
	}
	else if (sx < -1.0 + FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		xAngle = -HALF_PI;
		zAngle = 0.0;  // any angle works

		double fRmY = atan2(-m_data[0][1], m_data[0][0]);
		yAngle = zAngle - fRmY;
		return false;
	}

	xAngle = asin(sx);
	yAngle = atan2(m_data[0][2], m_data[2][2]);
	zAngle = atan2(m_data[1][0], m_data[1][1]);
	return true;
}

//-----------------------------------------------------------------------
bool CMatrix3x3::ToEulerAnglesYZX(double& yAngle, double& zAngle,
	double& xAngle) const
{
	// rot =  cy*cz           sx*sy-cx*cy*sz  cx*sy+cy*sx*sz
	//        sz              cx*cz          -cz*sx
	//       -cz*sy           cy*sx+cx*sy*sz  cx*cy-sx*sy*sz

	if (m_data[1][0] > 1.0 - FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		zAngle = HALF_PI;
		xAngle = 0.0;       // any angle works

		double fRpY = atan2(m_data[1][1], m_data[1][2]);
		yAngle = fRpY - xAngle;
		return false;
	}
	else if (m_data[1][0] < -1.0 + FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		zAngle = -HALF_PI;
		xAngle = 0.0;       // any angle works

		double fRmY = atan2(m_data[1][1], m_data[1][2]);
		yAngle = xAngle - fRmY;
		return false;
	}

	zAngle = asin(m_data[1][0]);
	yAngle = atan2(-m_data[2][0], m_data[0][0]);
	xAngle = atan2(-m_data[1][2], m_data[1][1]);
	return true;
}

//-----------------------------------------------------------------------
bool CMatrix3x3::ToEulerAnglesZXY(double& zAngle, double& xAngle,
	double& yAngle) const
{
	// rot =  cy*cz-sx*sy*sz -cx*sz           cz*sy+cy*sx*sz
	//        cz*sx*sy+cy*sz  cx*cz          -cy*cz*sx+sy*sz
	//       -cx*sy           sx              cx*cy

	if (m_data[2][1] > 1.0 - FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		xAngle = HALF_PI;
		yAngle = 0.0;        // any angle works

		double fRpY = atan2(m_data[0][2], m_data[0][0]);
		zAngle = fRpY - yAngle;
		return false;
	}
	else if (m_data[2][1] < -1.0 + FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		xAngle = -HALF_PI;
		yAngle = 0.0;        // any angle works

		double fRmY = atan2(m_data[0][2], m_data[0][0]);
		zAngle = yAngle - fRmY;
		return false;
	}

	xAngle = asin(m_data[2][1]);
	zAngle = atan2(-m_data[0][1], m_data[1][1]);
	yAngle = atan2(-m_data[2][0], m_data[2][2]);
	return true;
}
//-----------------------------------------------------------------------
bool CMatrix3x3::ToEulerAnglesZYX(double& zAngle, double& yAngle,
	double& xAngle) const
{
	// rot =  cy*cz           cz*sx*sy-cx*sz  cx*cz*sy+sx*sz
	//        cy*sz           cx*cz+sx*sy*sz -cz*sx+cx*sy*sz
	//       -sy              cy*sx           cx*cy

	double sy = -m_data[2][0];
	if (sy > 1.0 - FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		yAngle = HALF_PI;
		xAngle = 0.0;     // any angle works

		double fRpY = atan2(-m_data[0][1], m_data[0][2]);
		zAngle = fRpY - xAngle;
		return false;
	}
	else if (sy < -1.0 + FLT_EPSILON)
	{
		// WARNING.  Not a unique solution.
		yAngle = -HALF_PI;
		xAngle = 0.0;     // any angle works

		double fRmY = atan2(-m_data[0][1], m_data[0][2]);
		zAngle = xAngle - fRmY;
		return false;
	}

	yAngle = asin(sy);
	zAngle = atan2(m_data[1][0], m_data[0][0]);
	xAngle = atan2(m_data[2][1], m_data[2][2]);
	return true;
}

//-----------------------------------------------------------------------
void CMatrix3x3::FromEulerAnglesXYZ(const double& fYAngle, const double& fPAngle,
	const double& fRAngle)
{
	double fCos, fSin;

	fCos = cos(fYAngle);
	fSin = sin(fYAngle);
	CMatrix3x3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

	fCos = cos(fPAngle);
	fSin = sin(fPAngle);
	CMatrix3x3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

	fCos = cos(fRAngle);
	fSin = sin(fRAngle);
	CMatrix3x3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

	*this = kXMat*(kYMat*kZMat);
}
//-----------------------------------------------------------------------
void CMatrix3x3::FromEulerAnglesXZY(const double& fYAngle, const double& fPAngle,
	const double& fRAngle)
{
	double fCos, fSin;

	fCos = cos(fYAngle);
	fSin = sin(fYAngle);
	CMatrix3x3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

	fCos = cos(fPAngle);
	fSin = sin(fPAngle);
	CMatrix3x3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

	fCos = cos(fRAngle);
	fSin = sin(fRAngle);
	CMatrix3x3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

	*this = kXMat*(kZMat*kYMat);
}
//-----------------------------------------------------------------------
void CMatrix3x3::FromEulerAnglesYXZ(const double& fYAngle, const double& fPAngle,
	const double& fRAngle)
{
	double fCos, fSin;

	fCos = cos(fYAngle);
	fSin = sin(fYAngle);
	CMatrix3x3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

	fCos = cos(fPAngle);
	fSin = sin(fPAngle);
	CMatrix3x3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

	fCos = cos(fRAngle);
	fSin = sin(fRAngle);
	CMatrix3x3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

	*this = kYMat*(kXMat*kZMat);
}
//-----------------------------------------------------------------------
void CMatrix3x3::FromEulerAnglesYZX(const double& fYAngle, const double& fPAngle,
	const double& fRAngle)
{
	double fCos, fSin;

	fCos = cos(fYAngle);
	fSin = sin(fYAngle);
	CMatrix3x3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

	fCos = cos(fPAngle);
	fSin = sin(fPAngle);
	CMatrix3x3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

	fCos = cos(fRAngle);
	fSin = sin(fRAngle);
	CMatrix3x3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

	*this = kYMat*(kZMat*kXMat);
}
//-----------------------------------------------------------------------
void CMatrix3x3::FromEulerAnglesZXY(const double& fYAngle, const double& fPAngle,
	const double& fRAngle)
{
	double fCos, fSin;

	fCos = cos(fYAngle);
	fSin = sin(fYAngle);
	CMatrix3x3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

	fCos = cos(fPAngle);
	fSin = sin(fPAngle);
	CMatrix3x3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

	fCos = cos(fRAngle);
	fSin = sin(fRAngle);
	CMatrix3x3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

	*this = kZMat*(kXMat*kYMat);
}
//-----------------------------------------------------------------------
void CMatrix3x3::FromEulerAnglesZYX(const double& fYAngle, const double& fPAngle,
	const double& fRAngle)
{
	double fCos, fSin;

	fCos = cos(fYAngle);
	fSin = sin(fYAngle);
	CMatrix3x3 kZMat(fCos, -fSin, 0.0, fSin, fCos, 0.0, 0.0, 0.0, 1.0);

	fCos = cos(fPAngle);
	fSin = sin(fPAngle);
	CMatrix3x3 kYMat(fCos, 0.0, fSin, 0.0, 1.0, 0.0, -fSin, 0.0, fCos);

	fCos = cos(fRAngle);
	fSin = sin(fRAngle);
	CMatrix3x3 kXMat(1.0, 0.0, 0.0, 0.0, fCos, -fSin, 0.0, fSin, fCos);

	*this = kZMat*(kYMat*kXMat);
}

CVector3D operator *(const CMatrix3x3 &matrix, const CVector3D &vec)
{
#define MATRIX_MULT(i)  matrix[i][0] * vec[0] + matrix[i][1] * vec[1] + matrix[i][2] * vec[2]
	CVector3D result;
	result[0] = MATRIX_MULT(0);
	result[1] = MATRIX_MULT(1);
	result[2] = MATRIX_MULT(2);



	return result;
}

CMatrix3x3 CMatrix3x3::s_GetCovMatrix(const std::vector<CPoint3D> &vecPnts)
{
	int nPnts = vecPnts.size();
	if (nPnts < 3)
		return CMatrix3x3();

	CPoint3D mass;
	for (int i = 0; i < nPnts; i++)
	{
		mass += vecPnts[i];
	}
	mass /= nPnts;

	CMatrix3x3 matrix;
	for (int i = 0; i < nPnts; i++)
	{
		CPoint3D pnt = vecPnts[i] - mass;

		matrix[0][0] = pnt[0] * pnt[0];
		matrix[0][1] = pnt[0] * pnt[1];
		matrix[0][2] = pnt[0] * pnt[2];

		matrix[1][1] = pnt[1] * pnt[1];
		matrix[1][2] = pnt[1] * pnt[2];

		matrix[2][2] = pnt[2] * pnt[2];
	}

	matrix[1][0] = matrix[0][1];
	matrix[2][0] = matrix[0][2];

	matrix[2][1] = matrix[1][2];

	return matrix;
}

double CMatrix3x3::Determinant() const
{
	double det = m_data[0][0] * (m_data[1][1] * m_data[2][2] - m_data[1][2] * m_data[2][1])
		- m_data[1][0] * (m_data[0][1] * m_data[2][2] - m_data[0][2] * m_data[2][1])
		+ m_data[2][0] * (m_data[0][1] * m_data[1][2] - m_data[0][2] * m_data[1][1]);

	return det;
}

CMatrix3x3 CMatrix3x3::Inverse() const
{
	CMatrix3x3 result;

	double det = Determinant();
	if (fabs(det) < 1.0e-16)
	{
		return result;
	}

	det = 1.0 / det;

	//CMatrix3x3 companion_matrix;
	//std::cout << Determinant() <<"del"<< std::endl;
	result[0][0] = (m_data[1][1] * m_data[2][2] - m_data[1][2] * m_data[2][1]) * det;
	result[0][1] = (m_data[0][2] * m_data[2][1] - m_data[0][1] * m_data[2][2]) * det;
	result[0][2] = (m_data[0][1] * m_data[1][2] - m_data[0][2] * m_data[1][1]) * det;
	result[1][0] = (m_data[1][2] * m_data[2][0] - m_data[1][0] * m_data[2][2]) * det;
	result[1][1] = (m_data[0][0] * m_data[2][2] - m_data[0][2] * m_data[2][0]) * det;
	result[1][2] = (m_data[1][0] * m_data[0][2] - m_data[0][0] * m_data[1][2]) * det;
	result[2][0] = (m_data[1][0] * m_data[2][1] - m_data[1][1] * m_data[2][0]) * det;
	result[2][1] = (m_data[0][1] * m_data[2][0] - m_data[0][0] * m_data[2][1]) * det;
	result[2][2] = (m_data[0][0] * m_data[1][1] - m_data[1][0] * m_data[0][1]) * det;

	return result;
}

