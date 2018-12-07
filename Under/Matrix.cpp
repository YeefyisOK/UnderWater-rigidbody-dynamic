#include "math.h"
#include "Matrix.h"
//#include <math.h>

CMatrix::CMatrix(void)
{
	//memset((float *)m_data, sizeof(float) * 9, 0);
	for (int i = 0; i < 3; ++i)
	{
		for (int j = 0; j < 3; ++j)
		{
			m_data[i][j] = 0;
		}
	}
	m_data[0][0] = 1;
	m_data[1][1] = 1;
	m_data[2][2] = 1;
}


CMatrix::~CMatrix(void)
{
}

CVector3f CMatrix::operator *(const CVector3f &vector) const
{
	float res[3] = { 0, 0, 0 };

	res[0] = (m_data[0][0] + m_data[0][1] + m_data[0][2]) * vector[0];
	res[1] = (m_data[1][0] + m_data[1][1] + m_data[1][2]) * vector[1];
	res[2] = (m_data[2][0] + m_data[2][1] + m_data[2][2]) * vector[2];

	CVector3f v_result(res);
	return v_result;
}


CMatrix CMatrix::operator %(const CMatrix &matrix) const
{
#define MATRIX_MULT(i, j)  m_data[i][0] * matrix[0][j] + m_data[i][1] * matrix[1][j] + m_data[i][2] * matrix[2][j]
	CMatrix result;
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

CMatrix CMatrix::s_GetRotaionMatrix(float angle, const CPoint3D &axis)
{
	float x = axis[0];
	float y = axis[1];
	float z = axis[2];

#define  ZERO_TOL    1.0e-7
	float s = sin(angle);
	float c = cos(angle);

	CMatrix matrix;
	if (fabs(z) < ZERO_TOL)
	{
		if (fabs(y) < ZERO_TOL)
		{
			//if the axis is (0, 0, 0), assume to be identity matrix
			if (fabs(x) < ZERO_TOL)
				return matrix;

			//rotation axis is (1, 0, 0)
			matrix[0][0] = 1.0;
			matrix[1][1] = c;  matrix[2][2] = c;

			if (x > 0)
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

			if (y > 0)
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
	float mag = sqrt(x * x + y * y + z * z);
	mag = 1.0 / mag;
	x *= mag;
	y *= mag;
	z *= mag;

	float t = 1.0 - c;

	float tx = t * x;
	float ty = t * y;
	float tz = t * z;
	float sx = s * x;
	float sy = s * y;
	float sz = s * z;

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