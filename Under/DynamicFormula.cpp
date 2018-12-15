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
VectorXf DynamicFormula::tsfs2tf( Matrix3f Y) {
	Matrix3f Rt = R.transpose();
	cout << "RT" << Rt << endl;
	cout << "R-1" << R.inverse() << endl;
	Matrix3f zero = Matrix3f::Zero();
	//zero.setZero(3, 3);
	Matrix3f negRtY = zero - Rt * Y;
	cout << "negRtY" << negRtY << endl;
	MatrixXf trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt;
	trans.block(0, 3, 3, 3) = negRtY;
	trans.block(3, 0, 3, 3) = zero;
	trans.block(3, 3, 3, 3) = Rt;
	cout << "trans" << trans << endl;
	VectorXf tsfs(6);
	tsfs.block(0, 0, 3, 1) = ts;
	tsfs.block(3, 0, 3, 1) = fs;
	cout << "tsfs" << tsfs << endl;
	/*
	tsfs(0) = ts(0);
	tsfs(1) = ts(1);
	tsfs(2) = ts(2);
	tsfs(3) = fs(0);
	tsfs(4) = fs(1);
	tsfs(5) = fs(2);*/
	return trans * tsfs;
}
/*
Matrix3f DynamicFormula::computeR_() {
	Matrix3f daOmega = toDaOmegaOrY(w);
	return R * daOmega;
}*/
Vector3f DynamicFormula::computey_() {
	return R * v;
}
VectorXf DynamicFormula::computelp() {
	VectorXf wv(6);
	wv.block(0, 0, 3, 1) = w;
	wv.block(3, 0, 3, 1) = v;
	return K * wv;
}
VectorXf DynamicFormula::computelp_() {
	Matrix3f Y = toDaOmegaOrY(y);
	VectorXf tf = tsfs2tf(Y);
	cout << "tf" << tf << endl;
	Vector3f l = vec62Vec31(lp);
	Vector3f p = vec62Vec32(lp);
	cout << "L" << l << endl;
	cout << "P" << p << endl;
	Vector3f a = l.cross(w) + p.cross(v);
	cout << "A" << a << endl;
	Vector3f b = p.cross(w);
	cout << "B" << b << endl;
		VectorXf ab(6);
	ab.block(0, 0, 3, 1) = a;
	ab.block(3, 0, 3, 1) = b;
	return ab + tf;
}
Matrix3f DynamicFormula::computeNextR() {
	cout << "q:" << q.coeffs() << endl;
	float length = sqrt(w(0)*w(0) + w(1)*w(1) + w(2)*w(2));
	float q0 = cos(delta_t*length / 2);
	float genhao = sin(delta_t*length / 2);
	float q1, q2, q3;
	if(length != 0){
		q1 = genhao * w(0) / length;
		q2 = genhao * w(1) / length;
		q3 = genhao * w(2) / length;
	}
	else {
		q1 = q2 = q3 = 0;
	}
	Quaternionf tempq(q0, q1, q2, q3);
	delta_q = tempq;
	cout << "delta_q:" << delta_q.coeffs() << endl;
	temp_rotate = delta_q.toRotationMatrix();
	q = delta_q * q;//计算经过旋转后的orientation 
	q.normalize();
	return q.toRotationMatrix();
	/*
	cout << "q_old:" << q.coeffs() << endl;
	Quaternionf q_w(0, w(0)*delta_t, w(1)*delta_t, w(2)*delta_t);
	q_w *= q;
	q_w.w() = q_w.w()*0.5;
	q_w.x() = q_w.x()*0.5;
	q_w.y() = q_w.y()*0.5;
	q_w.z() = q_w.z()*0.5;
	Quaternionf delta_q=q_w;//以上是算delta_q
	cout << "delta_q;" << delta_q.coeffs() << endl;
	delta_q.normalize();
	temp_rotate = delta_q.toRotationMatrix();//delta_q转化为矩阵 display中与栈顶矩阵相乘（转4*4）
	cout <<"temp_rotate="<< temp_rotate << endl;
	q = delta_q * q ;//计算经过旋转后的orientation 
	q.normalize();
	R = q.toRotationMatrix();
	cout << "R:" << R << endl;*/
}
float* DynamicFormula::GetRotationData() {
	float data[16];
	data[0] = temp_rotate(0,0);
	data[1] = temp_rotate(1,0);
	data[2] = temp_rotate(2,0);
	data[3] = 0;

	data[4] = temp_rotate(0,1);
	data[5] = temp_rotate(1,1);
	data[6] = temp_rotate(2,1);
	data[7] = 0;

	data[8] = temp_rotate(0,2);
	data[9] = temp_rotate(1,2);
	data[10] = temp_rotate(2,2);
	data[11] = 0;

	data[12] = 0;
	data[13] = 0;
	data[14] = 0;
	data[15] = 1;

	return data;
}
Vector3f DynamicFormula::computeNexty( Vector3f y_) {
	temp_deltay = delta_t  * y_;
	cout << "(" << y(0) << "," << y(1) << "," << y(2) << ")" << endl;
	return y + temp_deltay;
}
VectorXf DynamicFormula::computeNextlp() {
	return lp + delta_t * lp_;
}
VectorXf DynamicFormula::computeNextwv() {
	MatrixXf Kneg = K.inverse();
	VectorXf res=Kneg * lp;
	//w = res.block(0, 0, 3, 1);
	//v = res.block(3, 0, 3, 1);
	return res;
}
Vector3f DynamicFormula::vec62Vec31(VectorXf wv) {//也可以用于wv
	return wv.block(0, 0, 3, 1);
}
Vector3f DynamicFormula::vec62Vec32(VectorXf wv) {
	return wv.block(3, 0, 3, 1);
}

void DynamicFormula::nextTime() {
	lp_ = computelp_();
	cout << "lp_:" << lp_ << endl;
	lp=computeNextlp();
	cout << "lp:" << lp << endl;
	VectorXf tempwv= computeNextwv();
	cout << "w:"<<w(0)<<" " << w(1) << " " << w(2) << endl;
	cout << "v:"<<v(0) << " " << v(1) << " " << v(2) << endl;
	Vector3f y_ = computey_();
	cout << "y_" << y_ << endl;
	y=computeNexty(y_);
	cout << "y:" << y << endl;
	R=computeNextR();
	w = tempwv.block(0, 0, 3, 1);
	v = tempwv.block(3, 0, 3, 1);
	cout << "w:" << w(0) << " " << w(1) << " " << w(2) << endl;
	cout << "v:" << v(0) << " " << v(1) << " " << v(2) << endl;
}
void DynamicFormula::set_tsfs(Vector3f ts,Vector3f fs) {
	this->ts = ts;
	this->fs = fs;
}