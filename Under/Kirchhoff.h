#pragma once
#include "PIC.h"
#include <vector>
#include"PICnew.h"
#include "Point3D.h"
#include "Vector6D.h"
#include <Eigen/Dense>
#include <Eigen/Cholesky>  
#include <Eigen/LU>  
#include <Eigen/QR>  
#include <Eigen/SVD>  
#include <math.h>
#include <vector>
using namespace Eigen;
class CKirchhoff
{
public:
	int numPoints;
	int numFaces;
	MatrixXd vertex ;//= MatrixXd::Random(3, 3)
	MatrixXd normal ;//= MatrixXd::Random(3, 3)
	MatrixXd face[3] ;
	double bodyMass;//computeJ得到了mass物体质量
	double fluidMass;//得到了流体质量
	double bodyDensity;
	double fluidDensity;
	double volume;
	Vector3d Cm;//质心，质量均匀分布
	CKirchhoff(PICnew *
		,double m_bodyDensity,double m_fluidDensity);
	MatrixXd angular_vector();
	double area();
	MatrixXd triangle_area();
	MatrixXd area_vector();
	MatrixXd face_center();
	MatrixXd face_normal();
	MatrixXd single_layer(MatrixXd a,MatrixXd b);
	MatrixXd motion_flux();
	MatrixXd one_point_quadrature();
	MatrixXd solid_angle(MatrixXd src);
	MatrixXd computeKF(double offset);
	//计算KB
	void Subexpressions(double &w0, double &w1, double &w2,
		double &f1, double &f2, double &f3, double &g0, double &g1, double &g2);
	Matrix3d computeJ();
	MatrixXd computeKB();
	MatrixXd computeK();//计算K=KF+KB
	VectorXd computetsfs();
};