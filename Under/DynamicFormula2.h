#pragma once
#include <Eigen/Dense>
#include <Eigen/Cholesky>  
#include <Eigen/LU>  
#include <Eigen/QR>  
#include <Eigen/SVD>  
//#include"Quaternion.h"
#include"Body.h"
#include <iostream>
using namespace std;
using namespace Eigen;
class DynamicFormula2
{
public:
	static VectorXd computetraction(vector<Body*> m_body);//计算每个面上的应力
	static Matrix3d computeKij(Vector3d x,Vector3d nx, Vector3d y,Vector3d ny);
	static Matrix3d computeHij(Vector3d x, Vector3d nx, Vector3d y, Vector3d ny);
	static double dirac(int a, int b);
	static double trianglearea(Vector3d p0,Vector3d p01,Vector3d p02) {
		double a0 = (p0 - p01).norm();
		double b0 = (p0 - p02).norm();
		double c0 = (p01 - p02).norm();
		double d0 = (a0 + b0 + c0) / 2;
		double S0 = sqrt(d0*(d0 - a0)*(d0 - b0)*(d0 - c0));
		return S0;
	}
/*
	MatrixXd K;
	double delta_t;*/
	static double mu;
	//平移旋转需要的两个参数
	//Quaternionf delta_q;
	//double theta;

};