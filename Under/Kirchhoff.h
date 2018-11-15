#pragma once
#include "PIC.h"
#include <vector>
#include "Point3D.h"
#include "Vector6D.h"
#include <Eigen/Cholesky>
#include <Eigen/Dense>
#include <math.h>
#include <vector>
class CKirchhoff
{
public:
	int numPoints;
	int numFaces;
	MatrixXd vertex;
	MatrixXd normal;
	MatrixXd face[3];

	CKirchhoff(PIC m_pic);
	~CKirchhoff(void);
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

	void Subexpressions(double &w0, double &w1, double &w2,
		double &f1, double &f2, double &f3, double &g0, double &g1, double &g2);
	MatrixXd computeKB(MatrixXd face[], int numFaces, int index[], double mass);
	CVector6D computeK();
};