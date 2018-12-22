#pragma once
#include "PIC.h"
#include <vector>
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
	MatrixXf vertex ;//= MatrixXf::Random(3, 3)
	MatrixXf normal ;//= MatrixXf::Random(3, 3)
	MatrixXf face[3] ;
	float bodyMass;//computeJ得到了mass
	float bodyDensity;
	float fluidDensity;
	float volume;

	CKirchhoff(PIC m_pic,float m_bodyDensity,float m_fluidDensity);
	MatrixXf angular_vector();
	float area();
	MatrixXf triangle_area();
	MatrixXf area_vector();
	MatrixXf face_center();
	MatrixXf face_normal();
	MatrixXf single_layer(MatrixXf a,MatrixXf b);
	MatrixXf motion_flux();
	MatrixXf one_point_quadrature();
	MatrixXf solid_angle(MatrixXf src);
	MatrixXf computeKF(float offset);
	//计算KB
	void Subexpressions(float &w0, float &w1, float &w2,
		float &f1, float &f2, float &f3, float &g0, float &g1, float &g2);
	Matrix3f comuputeJ();
	MatrixXf computeKB();
	MatrixXf computeK();//计算K=KF+KB
	VectorXf computetsfs();
};