#pragma once
#include <Eigen/Dense>
#include <Eigen/Cholesky>  
#include <Eigen/LU>  
#include <Eigen/QR>  
#include <Eigen/SVD>  
#include <vector>
#include "PICnew.h"
#include "PIC.h"
using namespace std;
using namespace Eigen;
struct onepointST
{
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	int id;
	Vector3d midpoint;//中点
	Vector3d normal;//法向
	double area;
};
class Body {
public:
	EIGEN_MAKE_ALIGNED_OPERATOR_NEW
	static int idnum;//把所有面片编号，用来赋值id
	Matrix4d g;
	VectorXd epsilon;
	VectorXd tf;
	//Vector3d f;
	MatrixXd K;
	int faceNum;
	double delta_t = 0.04;
	double bodyMass;//computeJ得到了mass
	double bodyDensity;//还没用到
	double fluidDensity;
	double volume;
	Vector3d masscenter;//世界坐标系的质心
	//VectorXd traction;
	PICnew *m_picnew;
	vector <onepointST> v_onepoint;//物体坐标变换到世界坐标
	Body(PICnew *m_picnew,Matrix3d R,Vector3d y);
	void Subexpressions(double &w0, double &w1, double &w2, double &f1, double &f2, double &f3, double &g0, double &g1, double &g2);
	Matrix3d computeJ();
	MatrixXd computeKB();

	float* GetRotAndTransData();
	void computetf(VectorXd traction);//通过所有面上的力计算合外力和合外力矩
	Matrix3d so3_ad(Vector3d omega);
	Matrix4d se3_cay(VectorXd tempepsilon);
	Matrix3d so3_cay(Vector3d tempw);
	MatrixXd se3_ad(VectorXd tempepsilon);
	VectorXd Unconstr_Dyn(VectorXd epsilon_now, VectorXd epsilon_last, Matrix4d &gk);
	VectorXd se3_DEP(VectorXd epsilon_now, VectorXd epsilon_last, Matrix4d &gk);
	MatrixXd se3_Ctln(VectorXd tempepsilon);
	void nextTime();
	
};