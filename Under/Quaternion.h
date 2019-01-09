#pragma once
#include "Matrix.h"
#include <math.h>

class CQuaternion
{
public:
	CQuaternion(double _w=0,double _x=1,double _y=0,double _z=0)
	{
		w = _w;
		x = _x;
		y = _y;
		z = _z;
	}
	CMatrix ToMatrix()
	{
		CMatrix mat;
		double xx2 = 2 * x*x, yy2 = 2 * y*y, zz2 = 2 * z*z,
			xy2 = 2 * x*y, xz2 = 2 * x*z, yz2 = 2 * y*z,
			xw2 = 2 * x*w, yw2 = 2 * y*w, zw2 = 2 * z*w;
		mat[0][0] = 1 - yy2 - zz2;
		mat[1][0] = xy2 + zw2;
		mat[2][0] = xz2 - yw2;
		mat[0][1] = xy2 - zw2;
		mat[1][1] = 1 - xx2 - zz2;
		mat[2][1] = yz2 + xw2;
		mat[0][2] = xz2 + yw2;
		mat[1][2] = yz2 - xw2;
		mat[2][2] = 1 - xx2 - yy2;
		return mat;
	}

	void FromRotationMatrix(const CMatrix& mat)
	{
		double s;
		double const tr = mat[0][0] + mat[1][1] + mat[2][2] + 1.f;

		// check the diagonal
		if (tr > 1.f)
		{
			s = sqrt(tr);
			w = s * 0.5f;
			s = 0.5f / s;
			x = (mat[1][2] - mat[2][1]) * s;
			y = (mat[2][0] - mat[0][2]) * s;
			z = (mat[0][1] - mat[1][0]) * s;
		}
		else
		{
			int maxi = 0;
			double maxdiag = mat[0][0];
			for (int i = 1; i < 3; ++i)
			{
				if (mat[i][i] > maxdiag)
				{
					maxi = i;
					maxdiag = mat[i][i];
				}
			}

			switch (maxi)
			{
			case 0:
				s = sqrt((mat[0][0] - (mat[1][1] + mat[2][2])) + 1);

				x = s * 0.5f;

				if (s>0.05f)
				{
					s = 0.5f / s;
				}
				w = (mat[1][2] - mat[2][1]) * s;
				y = (mat[1][0] + mat[0][1]) * s;
				z = (mat[2][0] + mat[0][2]) * s;
				break;

			case 1:
				s = sqrt((mat[1][1] - (mat[2][2] + mat[0][0])) + 1);
				y = s * 0.5f;

				if (s>0.05f)
				{
					s = 0.5f / s;
				}

				w = (mat[2][0] - mat[0][2]) * s;
				z = (mat[2][1] + mat[1][2]) * s;
				x = (mat[0][1] + mat[1][0]) * s;
				break;

			case 2:
			default:
				s = sqrt((mat[2][2] - (mat[0][0] + mat[1][1])) + 1);
				z = s * 0.5f;
				if (s>0.05f)
				{
					s = 0.5f / s;
				}

				w = (mat[0][1] - mat[1][0]) * s;
				x = (mat[0][2] + mat[2][0]) * s;
				y = (mat[1][2] + mat[2][1]) * s;
				break;
			}
		}
	}
	static CQuaternion FromRotation(double _x,double _y,double _z,double angleInDegree)
	{
		// 向量的单位化
		double length = sqrt(_x * _x + _y * _y + _z * _z);
		_x /= length;
		_y /= length;
		_z /= length;
		double alpha = angleInDegree / 180 * 3.1415926;// 已转换弧度制
		return CQuaternion(
			sin(alpha / 2) * _x,
			sin(alpha / 2) * _y,
			sin(alpha / 2) * _z,
			cos(alpha / 2));
	}
	void setomega(double new_w){
		w = new_w;
	}
	double get_w(){
		return w;
	}
	double get_x(){
		return x;
	}
	double get_y(){
		return y;
	}
	double get_z(){
		return z;
	}
private:
	double  w, x, y, z;
};