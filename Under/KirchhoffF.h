#pragma once
#include <Eigen/Dense>
#include "PIC.h"
#include<vector>
class CKirchhoffF
{
public:
	int numPoints;
	int numFaces;
	MatrixXd vertex;
	MatrixXd normal;
	MatrixXd face[3];

	CKirchhoffF(PIC m_pic);
	~CKirchhoffF(void);
	MatrixXd angular_vector();
	MatrixXd area();
	MatrixXd triangle_area();
	MatrixXd area_vector();
	double surf_area();
	MatrixXd face_center();
	MatrixXd face_normal();
	MatrixXd single_layer();
	MatrixXd double_layer();
	MatrixXd motion_flux();
	MatrixXd one_point_quadrature();
	MatrixXd solid_angle(MatrixXd src);
	MatrixXd compute(double offset);
};