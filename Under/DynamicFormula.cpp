#include"DynamicFormula.h"


DynamicFormula::DynamicFormula(Vector3d omega, Vector3d velocity, Matrix3d R,
	Vector3d y, Vector3d ts, Vector3d fs, MatrixXd K, double delta_t)
{
	this->w = omega;
	this->v = velocity;
	this->R = R; //R如何初始化
	this->y = y;
	this->ts = ts;
	this->fs = fs;
	this->K = K;
	this->delta_t = delta_t;
}

DynamicFormula::~DynamicFormula()
{
}
Matrix3d DynamicFormula::toDaOmegaOrY(Vector3d omega) {
	Matrix3d res;
	res(0, 0) = 0;
	res(0, 1) = -omega(2);
	res(0, 2) = omega(1);
	res(1, 0) = omega(2);
	res(1, 1) = 0;
	res(1, 2) = -omega(0);
	res(2, 0) = -omega(1);
	res(2, 1) = omega(0);
	res(2, 2) = 0;
	return res;
}
VectorXd DynamicFormula::tsfs2tf(Matrix3d R, Matrix3d Y, Vector3d ts, Vector3d fs){
	Matrix3d Rt = R.transpose();
	Matrix3d zero=Matrix3d::Zero();
	//zero.setZero(3, 3);
	Matrix3d negRtY =  zero- Rt * Y;//????
	MatrixXd trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt.block(0,0,3,3);
	trans.block(0, 3, 3, 3) = negRtY.block(0, 0, 3, 3);
	trans.block(3, 0, 3, 3) = zero.block(0, 0, 3, 3);
	trans.block(3, 3, 3, 3) = Rt.block(0, 0, 3, 3);
	VectorXd tsfs(6);
	tsfs.block(0, 0, 3, 1) = ts;
	tsfs.block(3, 0, 3, 1) = fs;
	/*
	tsfs(0) = ts(0);
	tsfs(1) = ts(1);
	tsfs(2) = ts(2);
	tsfs(3) = fs(0);
	tsfs(4) = fs(1);
	tsfs(5) = fs(2);*/	
	return trans * tsfs;
}
Matrix3d DynamicFormula::computeR_() {
	Matrix3d daOmega = toDaOmegaOrY(w);
	return R * daOmega;
}
Vector3d DynamicFormula::computey_() {
	return R * v;
}
VectorXd DynamicFormula::computelp() {
	VectorXd wv(6);
	wv.block(0, 0, 3, 1) = w;
	wv.block(3, 0, 3, 1) = v;
	return K * wv;
}
VectorXd DynamicFormula::computelp_(VectorXd lp) {
	Matrix3d Y = toDaOmegaOrY(y);
	VectorXd tf= tsfs2tf(R, Y,ts, fs);
	Vector3d l = vec62Vec31(lp);
	Vector3d p = vec62Vec32(lp);
	Vector3d a = l.cross(w) + p.cross(v);
	Vector3d b = p.cross(w);
	VectorXd ab(6);
	ab.block(0, 0, 3, 1) = a;
	ab.block(3, 0, 3, 1) = b;
	return ab + tf;
}
void DynamicFormula::computeNextR(Matrix3d R_) {
	q = R_;//直接可以赋值
	Quaterniond newq(q.w()* delta_t, q.x(), q.y(), q.z());
	Matrix3d temp = newq.matrix();
	R= R * temp;
}
void DynamicFormula::computeNexty( Vector3d y_) {
	temp_deltay = delta_t * y_;
	y= y + temp_deltay;
	cout << "(" << y(0) << "," << y(1) << "," << y(2) << ")" << endl;
}
VectorXd DynamicFormula::computeNextlp(VectorXd lp, VectorXd lp_) {
	return lp + delta_t * lp_;
}
void DynamicFormula::computeNextwv(VectorXd lp) {
	MatrixXd Kneg = K.inverse();
	VectorXd res=Kneg * lp;
	w = res.block(0, 0, 3, 1);
	v = res.block(3, 0, 3, 1);
}
Vector3d DynamicFormula::vec62Vec31(VectorXd wv) {//也可以用于wv
	return wv.block(0, 0, 3, 1);
}
Vector3d DynamicFormula::vec62Vec32(VectorXd wv) {
	return wv.block(3, 0, 3, 1);
}

void DynamicFormula::nextTime() {
	Matrix3d R_ = computeR_();
	Vector3d y_ = computey_();
	VectorXd lp = computelp();
	VectorXd lp_ = computelp_(lp);
	computeNextR(R_);
	computeNexty(y_);
	lp=computeNextlp(lp,lp_);
	computeNextwv(lp);
}
void DynamicFormula::set_tsfs(Vector3d ts,Vector3d fs) {
	this->ts = ts;
	this->fs = fs;
}