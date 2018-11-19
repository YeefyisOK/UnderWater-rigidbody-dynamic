#pragma once
#include <Eigen/Dense>
#include <Eigen/Cholesky>  
#include <Eigen/LU>  
#include <Eigen/QR>  
#include <Eigen/SVD>  
#include"Quaternion.h"
using namespace Eigen;
class DynamicFormula
{
public:
	DynamicFormula(Vector3d omega, Vector3d velocity, Matrix3d R,
		Vector3d y, Vector3d ts, Vector3d fs,MatrixXd K,double delta_t);
	~DynamicFormula();
	void setK(MatrixXd newK) {
		K = newK;
	}
	Matrix3d toDaOmegaOrY(Vector3d omega);
	VectorXd tsfs2tf(Matrix3d R, Matrix3d Y, Vector3d ts, Vector3d fs);
	Matrix3d computeR_();
	Matrix3d computey_();
	VectorXd computelp();

	VectorXd computelp_(VectorXd lp);
	void computeNextR(Matrix3d R_);
	void computeNexty(Vector3d y_);
	VectorXd computeNextlp(VectorXd lp, VectorXd lp_);
	void computeNextwv(VectorXd lp);
	Vector3d vec62Vec31(VectorXd ab);
	Vector3d vec62Vec32(VectorXd ab);
	void nextTime();
	void set_tsfs(Vector3d ts,Vector3d fs);

//private:
	//初始时刻的 6+1个量
	Vector3d w;
	Vector3d v;
	Matrix3d R; //R初始化一个正交矩阵
	Vector3d y;
	Vector3d ts;
	Vector3d fs;

	MatrixXd K;
	double delta_t;
	//平移旋转需要的两个参数
	Quaterniond q;
	Vector3d temp_deltay;

};

DynamicFormula::DynamicFormula(Vector3d omega, Vector3d velocity,Matrix3d R,
	Vector3d y, Vector3d ts, Vector3d fs,MatrixXd K,double delta_t)
{
	this->w= omega;
	this->v= velocity;
	this->R	=R; //R如何初始化
	this->y=y;
	this->ts=ts;
	this->fs = fs;
	this->K = K;
	this->delta_t = delta_t;
}

DynamicFormula::~DynamicFormula()
{
}