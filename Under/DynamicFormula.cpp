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
	Matrix3f zero = Matrix3f::Zero();
	//zero.setZero(3, 3);
	Matrix3f negRtY = zero - Rt * Y;
	MatrixXf trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt;
	trans.block(0, 3, 3, 3) = negRtY;
	trans.block(3, 0, 3, 3) = zero;
	trans.block(3, 3, 3, 3) = Rt;
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
	Vector3f l = vec62Vec31(lp);
	Vector3f p = vec62Vec32(lp);
	Vector3f a = l.cross(w) + p.cross(v);
	Vector3f b = p.cross(w);
	VectorXf ab(6);
	ab.block(0, 0, 3, 1) = a;
	ab.block(3, 0, 3, 1) = b;
	return ab + tf;
}
Matrix3f DynamicFormula::s_GetRotaionMatrix(float angle, const Vector3f &axis)
{
	float x = axis(0);
	float y = axis(1);
	float z = axis(2);

#define  ZERO_TOL    1.0e-7
	float s = sin(angle);
	float c = cos(angle);

	Matrix3f matrix = Matrix3f::Identity();
	if (fabs(z) < ZERO_TOL)
	{
		if (fabs(y) < ZERO_TOL)
		{
			//if the axis is (0, 0, 0), assume to be identity matrix
			if (fabs(x) < ZERO_TOL)
				return matrix;

			//rotation axis is (1, 0, 0)
			matrix(0,0)= 1.0;
			matrix(1,1)= c;  matrix(2,2) = c;

			if (x > 0)
			{
				matrix(1,2)= -s; matrix(2,1) = s;
			}
			else
			{
				matrix(1,2) = s; matrix(2,1)= -s;
			}
			return matrix;
		}
		else if (fabs(x) < ZERO_TOL)
		{
			//rotation axis is (0, 1, 0)
			matrix(1,1)= 1.0;
			matrix(0,0)= c;  matrix(2,2)= c;
			if (y > 0)
			{
				matrix(0,2) = s;  matrix(2,0)= -s;
			}
			else
			{
				matrix(0,2) = -s;  matrix(2,0)= s;
			}
			return matrix;
		}
	}
	else if (fabs(y) < ZERO_TOL)
	{
		if (fabs(x) < ZERO_TOL)
		{
			//rotation axis is (0, 0, 1)
			matrix(2,2) = 1.0;
			matrix(0,0) = c;  matrix(1,1)= c;

			if (z > 0)
			{
				matrix(0,1)= -s; matrix(1,0) = s;
			}
			else
			{
				matrix(0,1)= s; matrix(1,0)= -s;
			}
			return matrix;
		}
	}
	//common case
	//normalize the rotation axis
	float mag = sqrt(x * x + y * y + z * z);
	mag = 1.0 / mag;
	x *= mag;
	y *= mag;
	z *= mag;
	float t = 1.0 - c;

	float tx = t * x;
	float ty = t * y;
	float tz = t * z;
	float sx = s * x;
	float sy = s * y;
	float sz = s * z;
	//-----------------------------------------------------------
	//		| t*x*x + c		t*x*y - s*z		t*x*z + s*y |
	//		|											|
	//	R = | t*x*y + s*z	t*y*y + c		t*y*z - s*x |
	//		|											|
	//		| t*x*z - s*y	t*y*z + s*x		t*z*z + c	|
	//
	// where c = cos(theta), s = sin(theta), t = 1 - c and(x, y, z) is a unit
	// vector on the axis of rotation.
	//-----------------------------------------------------------

	// row one
	matrix(0,0)= tx * x + c;
	matrix(0,1) = tx * y - sz;
	matrix(0,2)= tx * z + sy;

	// row two
	matrix(1,0) = matrix(0,1) + sz + sz;	// tx * y + sz
	matrix(1,1) = ty * y + c;
	matrix(1,2) = ty * z - sx;

	// row three
	matrix(2,0)= matrix(0,2) - sy - sy;	// tx * z - sy
	matrix(2,1) = matrix(1,2) + sx + sx;	// ty * z + sx
	matrix(2,2) = tz * z + c;
	return matrix;
}
Matrix3f DynamicFormula::computeNextR() {
	/*
	Matrix3f daomega = toDaOmegaOrY(w);
	daomega(0, 0) = daomega(0, 0) * delta_t+1;
	daomega(0, 1) *= delta_t;
	daomega(0, 2) *= delta_t;
	daomega(1, 0) *= delta_t;
	daomega(1, 1) = daomega(1, 1) * delta_t + 1;
	daomega(1, 2) *= delta_t;
	daomega(2, 0) *= delta_t;
	daomega(2, 1) *= delta_t;
	daomega(2, 2) = daomega(2, 2) * delta_t + 1;
	return daomega * R;*/
	/*
	Vector3f Rw = R * w;
	cout << "RW" << Rw << endl;
	float length = sqrt(Rw(0)*Rw(0) + Rw(1)*Rw(1) + Rw(2)*Rw(2));
	float angle = delta_t * length;
	Matrix3f delta_matrix = s_GetRotaionMatrix(angle, Rw);
	cout << "delta_matrix"<<delta_matrix << endl;
	return delta_matrix * R;*/
	/*
	cout << "q:" << q.coeffs() << endl;
	Vector3f Rw =R* w;
	cout << "RW" << Rw << endl;
	float length = sqrt(Rw(0)*Rw(0) + Rw(1)*Rw(1) + Rw(2)*Rw(2));
	theta = delta_t * length;
	float q0 = cos(delta_t*length / 2);
	float genhao = sin(delta_t*length / 2);
	float q1, q2, q3;
	if(length != 0){
		q1 = genhao * Rw(0) / length;
		q2 = genhao * Rw(1) / length;
		q3 = genhao * Rw(2) / length;
	}
	else {
		q1 = q2 = q3 = 0;
	}
	Quaternionf tempq(q0, q1, q2, q3);
	Matrix3f temR=tempq.toRotationMatrix();
	//q = tempq * q;//计算经过旋转后的orientation 
	//q.normalize();
	//return q.toRotationMatrix();
	return temR * R;*/
	/*
	cout << "q_old:" << q.coeffs() << endl;
	Vector3f Rw = R * w;
	Quaternionf q_w(0, Rw(0)*delta_t, Rw(1)*delta_t, Rw(2)*delta_t);
	q_w = q_w*q;
	q_w.w() = q_w.w()*0.5;
	q_w.x() = q_w.x()*0.5;
	q_w.y() = q_w.y()*0.5;
	q_w.z() = q_w.z()*0.5;
	cout << "q_w:" << q_w.coeffs() << endl;//计算经过旋转后的orientation 
	q = q_w * q;
	q.normalize();
	return q.toRotationMatrix();*/
	Vector3f Rw = R *w;// 
	float length = 1;//sqrt(Rw(0)*Rw(0) + Rw(1)*Rw(1) + Rw(2)*Rw(2))
	Quaternionf q_w(0, Rw(0)/ length *delta_t, Rw(1)/ length *delta_t, Rw(2)/ length *delta_t);
	q_w = q_w * q;
	q_w.w() = q_w.w()*0.5;
	q_w.x() = q_w.x()*0.5;
	q_w.y() = q_w.y()*0.5;
	q_w.z() = q_w.z()*0.5;
	cout << "q_w:" << q_w.coeffs() << endl;//计算经过旋转后的orientation 
	q.w() = q_w.w() + q.w();
	q.x() = q_w.x() + q.x();
	q.y() = q_w.y() + q.y();
	q.z() = q_w.z() + q.z();
	q.normalize();
	return q.toRotationMatrix();
}
float* DynamicFormula::GetRotAndTransData() {
	static float data[16];//!!!!
	data[0] = R(0,0);
	data[1] = R(1,0);
	data[2] = R(2,0);
	data[3] = 0;

	data[4] = R(0,1);
	data[5] = R(1,1);
	data[6] = R(2,1);
	data[7] = 0;

	data[8] = R(0,2);
	data[9] = R(1,2);
	data[10] = R(2,2);
	data[11] = 0;

	data[12] = 0;
	data[13] = 0;
	data[14] = 0;
	data[15] = 1;
	return data;
}
Vector3f DynamicFormula::computeNexty( Vector3f y_) {
	Vector3f temp_deltay = delta_t  * y_;
//	cout << "temp_deltay" << temp_deltay << endl;
//	cout << "(" << y(0) << "," << y(1) << "," << y(2) << ")" << endl;
	return y + temp_deltay;
}
VectorXf DynamicFormula::computeNextlp() {
	return lp + delta_t * lp_;
}
VectorXf DynamicFormula::computeNextwv() {
	MatrixXf Kinv = K.inverse();
	cout << "Kinv" << Kinv << endl;
	VectorXf res=Kinv * lp;
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
	cout << "oldlp" << lp << endl;
	cout << "lp_:" << lp_ << endl;
	lp=computeNextlp();//lp
	cout << "lp:" << lp << endl;
	VectorXf tempwv= computeNextwv();
	cout << "tempwv" << tempwv << endl;
	cout << "w:"<<w(0)<<" " << w(1) << " " << w(2) << endl;
	cout << "v:"<<v(0) << " " << v(1) << " " << v(2) << endl;
	Vector3f y_ = computey_();
	cout << "y_" << y_ << endl;
	y=computeNexty(y_);//y
	cout << "y:" << y << endl;
	R=computeNextR();//R
	cout << "R:" << R << endl;
	w = tempwv.block(0, 0, 3, 1);//w
	v = tempwv.block(3, 0, 3, 1);//v
	cout << "w:" << w(0) << " " << w(1) << " " << w(2) << endl;
	cout << "v:" << v(0) << " " << v(1) << " " << v(2) << endl;
}
void DynamicFormula::set_tsfs(Vector3f ts,Vector3f fs) {
	this->ts = ts;
	this->fs = fs;
}