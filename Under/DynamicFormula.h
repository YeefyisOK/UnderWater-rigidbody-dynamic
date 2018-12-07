#pragma once
#include <Eigen/Dense>
#include <Eigen/Cholesky>  
#include <Eigen/LU>  
#include <Eigen/QR>  
#include <Eigen/SVD>  
//#include"Quaternion.h"
#include <iostream>
using namespace std;
using namespace Eigen;
class DynamicFormula
{
public:
	DynamicFormula(Vector3f omega, Vector3f velocity, Matrix3f R,
		Vector3f y, Vector3f ts, Vector3f fs,MatrixXf K,float delta_t);
	~DynamicFormula();
	void setK(MatrixXf newK) {
		K = newK;
	}
	void setR(Matrix3f newR) {
		R = newR;
	}
	Matrix3f toDaOmegaOrY(Vector3f omega);
	VectorXf tsfs2tf(Matrix3f R, Matrix3f Y);
	Matrix3f computeR_();
	Vector3f computey_();
	VectorXf computelp();

	VectorXf computelp_(VectorXf lp);
	void computeNextR(Matrix3f R_);
	void computeNexty(Vector3f y_);
	VectorXf computeNextlp(VectorXf lp, VectorXf lp_);
	void computeNextwv(VectorXf lp);
	Vector3f vec62Vec31(VectorXf ab);
	Vector3f vec62Vec32(VectorXf ab);
	void nextTime();
	void set_tsfs(Vector3f ts,Vector3f fs);

//private:
	//初始时刻的 6+1个量
	Vector3f w;
	Vector3f v;
	Matrix3f R; //R初始化一个正交矩阵
	Vector3f y;
	Vector3f ts;
	Vector3f fs;

	MatrixXf K;
	float delta_t;
	//平移旋转需要的两个参数
	Quaternionf q;
	Vector3f temp_deltay;

};