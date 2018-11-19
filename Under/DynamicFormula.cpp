#include"DynamicFormula.h"

Matrix3d DynamicFormula::toDaOmegaOrY(Vector3d omega) {
	Matrix3d res;
	res(0, 0) = 0;
	res(0, 1) = -omega[2];
	res(0, 2) = omega[1];
	res(1, 0) = omega[2];
	res(1, 1) = 0;
	res(1, 2) = -omega[0];
	res(2, 0) = -omega[1];
	res(2, 1) = omega[0];
	res(2, 2) = 0;
	return res;
}
VectorXd DynamicFormula::tsfs2tf(Matrix3d R, Matrix3d Y, Vector3d ts, Vector3d fs){
	Matrix3d Rt = R.transpose();
	Matrix3d zero=Matrix3d::Zero();
	//zero.setZero(3, 3);
	Matrix3d negRtY = zero - Rt * Y;
	MatrixXd trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt;
	trans.block(0, 3, 3, 3) = negRtY;
	trans.block(3, 0, 3, 3) = zero;
	trans.block(3, 3, 3, 3) = Rt;
	VectorXd tsfs(6);
	tsfs.block(0, 0, 1, 3) = ts;
	tsfs.block(0, 3, 1, 3) = fs;
	VectorXd tf = trans * tsfs;
	return tf;
}
Matrix3d DynamicFormula::computeR_() {
	Matrix3d daOmega = toDaOmegaOrY(w);
	return R * daOmega;
}
Matrix3d DynamicFormula::computey_() {
	return R * v;
}
VectorXd DynamicFormula::computelp() {
	VectorXd wv(6);
	wv.block(0, 0, 1, 3) = w;
	wv.block(0, 3, 1, 3) = v;
	return K * wv;
}
VectorXd DynamicFormula::computelp_(VectorXd lp) {
	Matrix3d Y = toDaOmegaOrY(y);
	VectorXd tf=tsfs2tf(R, Y,ts, fs);
	Vector3d l = vec62Vec31(lp);
	Vector3d p = vec62Vec32(lp);
	Vector3d a = l.cross(w) + p.cross(v);
	Vector3d b = p.cross(w);
	VectorXd ab(6);
	ab.block(0, 0, 1, 3) = a;
	ab.block(0, 3, 1, 3) = b;
	return ab + tf;
}
void DynamicFormula::computeNextR(Matrix3d R_) {
	q = R_;//直接可以赋值
	Quaterniond newq(q.w()* delta_t, q.x(), q.y(), q.z());
	q = newq;
	Matrix3d temp = q.toRotationMatrix();
	R= R + temp;
}
void DynamicFormula::computeNexty( Vector3d y_) {
	temp_deltay = delta_t * y_;
	y= y + temp_deltay;
}
VectorXd DynamicFormula::computeNextlp(VectorXd lp, VectorXd lp_) {
	return lp + delta_t * lp_;
}
void DynamicFormula::computeNextwv(VectorXd lp) {
	Matrix3d Kneg = K.inverse();
	VectorXd res=Kneg * lp;
	w = res.block(0, 0, 1, 3);
	v = res.block(0, 3, 1, 3);
}
Vector3d DynamicFormula::vec62Vec31(VectorXd wv) {//也可以用于wv
	return wv.block(0, 0, 1, 3);
}
Vector3d DynamicFormula::vec62Vec32(VectorXd wv) {
	return wv.block(0, 3, 1, 3);
}

void DynamicFormula::nextTime() {
	Matrix3d R_ =computeR_();
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