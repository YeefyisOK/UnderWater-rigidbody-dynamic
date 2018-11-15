#pragma once
#include <Eigen/Dense>
#include "PIC.h"
#include<vector>
using namespace Eigen;
using namespace std;
class CKirchhoffB
{
public:
	int numPoints;
	int numFaces;
	MatrixXd vertex;
	MatrixXd normal;
	MatrixXd face[3];

	CKirchhoffB(PIC m_pic);
	MatrixXd compute(double offset);
	double surf_area();
	MatrixXd face_center();
	MatrixXd motion_flux();
	MatrixXd division(MatrixXd MF, MatrixXd M);
	MatrixXd one_point_quadrature();
	MatrixXd angular_vector();
	MatrixXd area();
	MatrixXd triangle_area();
	MatrixXd area_vector();
	MatrixXd solid_angle();

	~CKirchhoffB(void);
};