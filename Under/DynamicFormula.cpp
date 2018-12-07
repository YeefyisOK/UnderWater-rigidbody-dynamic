#include"DynamicFormula.h"
DynamicFormula::DynamicFormula(Vector3f omega, Vector3f velocity, Matrix3f R,
	Vector3f y, Vector3f ts, Vector3f fs, MatrixXf K, float delta_t)
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
Matrix3f DynamicFormula::toDaOmegaOrY(Vector3f omega) {
	Matrix3f res;
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
VectorXf DynamicFormula::tsfs2tf(Matrix3f R, Matrix3f Y){
	Matrix3f Rt = R.transpose();
	Matrix3f zero=Matrix3f::Zero();
	//zero.setZero(3, 3);
	Matrix3f negRtY =  zero- Rt * Y;
	MatrixXf trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt.block(0, 0, 3, 3);
	trans.block(0, 3, 3, 3) = negRtY.block(0, 0, 3, 3);
	trans.block(3, 0, 3, 3) = zero.block(0, 0, 3, 3);
	trans.block(3, 3, 3, 3) = Rt.block(0, 0, 3, 3);
	VectorXf tsfs(6);
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
Matrix3f DynamicFormula::computeR_() {
	Matrix3f daOmega = toDaOmegaOrY(w);
	return R * daOmega;
}
Vector3f DynamicFormula::computey_() {
	return R * v;
}
VectorXf DynamicFormula::computelp() {
	VectorXf wv(6);
	wv.block(0, 0, 3, 1) = w;
	wv.block(3, 0, 3, 1) = v;
	return K * wv;
}
VectorXf DynamicFormula::computelp_(VectorXf lp) {
	Matrix3f Y = toDaOmegaOrY(y);
	VectorXf tf= tsfs2tf(R, Y);
	Vector3f l = vec62Vec31(lp);
	Vector3f p = vec62Vec32(lp);
	Vector3f a = l.cross(w) + p.cross(v);
	Vector3f b = p.cross(w);
	VectorXf ab(6);
	ab.block(0, 0, 3, 1) = a;
	ab.block(3, 0, 3, 1) = b;
	return ab + tf;
}
void DynamicFormula::computeNextR(Matrix3f R_) {
//	cout << "R_=" << R_ << endl;
//	cout << "R=" << R << endl;
	Quaternionf q_old(R);//可以构造函数 也可以直接赋值
	q_old.normalize();
//	cout << "q_old:" << q_old.coeffs() << endl;
	Quaternionf q_w(0, w(0), w(1), w(2));
	float length =sqrt( w(0)*w(0) + w(1)*w(1) + w(2)*w(2));
	q_w.normalize();
	Quaternionf delta_q((q_w.w()*delta_t + 1)*length, q_w.x()*delta_t, q_w.y()*delta_t, q_w.z()*delta_t);
//	cout << "delta_q:" << delta_q.coeffs() << endl;
	q = delta_q;
//	cout << "delta_q:" << delta_q.coeffs() << endl;
	Quaternionf q_new = delta_q * q_old;
//	cout << "q_new:" << q_new.coeffs() << endl;
	q_new.normalize();
	R = q_new.toRotationMatrix();
}
void DynamicFormula::computeNexty( Vector3f y_) {
	temp_deltay = delta_t * y_;
	y= y + temp_deltay;
	cout << "(" << y(0) << "," << y(1) << "," << y(2) << ")" << endl;
}
VectorXf DynamicFormula::computeNextlp(VectorXf lp, VectorXf lp_) {
	return lp + delta_t * lp_;
}
void DynamicFormula::computeNextwv(VectorXf lp) {
	MatrixXf Kneg = K.inverse();
	VectorXf res=Kneg * lp;
	w = res.block(0, 0, 3, 1);
	v = res.block(3, 0, 3, 1);
}
Vector3f DynamicFormula::vec62Vec31(VectorXf wv) {//也可以用于wv
	return wv.block(0, 0, 3, 1);
}
Vector3f DynamicFormula::vec62Vec32(VectorXf wv) {
	return wv.block(3, 0, 3, 1);
}

void DynamicFormula::nextTime() {
	Matrix3f R_ = computeR_();
	Vector3f y_ = computey_();
	VectorXf lp = computelp();
	VectorXf lp_ = computelp_(lp);
	computeNextR(R_);
	computeNexty(y_);
	lp=computeNextlp(lp,lp_);
	computeNextwv(lp);
}
void DynamicFormula::set_tsfs(Vector3f ts,Vector3f fs) {
	this->ts = ts;
	this->fs = fs;
}