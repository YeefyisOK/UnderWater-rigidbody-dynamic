#include"DynamicFormula.h"
DynamicFormula::DynamicFormula(Vector3d omega, Vector3d velocity, Matrix3d R,
	Vector3d y,double delta_t)
{
	VectorXd temp(6);
	temp.block(0, 0, 3, 1) = omega;
	temp.block(3, 0, 3, 1) = velocity;
	epsilon = temp;
	this->w = omega;
	this->v = velocity;
	this->R = R; //R如何初始化
	this->y = y;
	this->delta_t = delta_t;
	g.block(0, 0, 3, 3) = R;
	g.block(0, 3, 3, 1) = y;
	g(3, 0) = 0;
	g(3, 1) = 0;
	g(3, 2) = 0;
	g(3, 3) = 1;
}

DynamicFormula::~DynamicFormula()
{
}
Matrix3d DynamicFormula::so3_ad(Vector3d omega) {
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
VectorXd DynamicFormula::tsfs2tf( Matrix3d Y) {
	Matrix3d Rt = R.transpose();
	Matrix3d zero = Matrix3d::Zero();
	//zero.setZero(3, 3);
	Matrix3d negRtY = zero - Rt * Y;
	MatrixXd trans(6, 6);//矩阵分块赋值
	trans.block(0, 0, 3, 3) = Rt;
	trans.block(0, 3, 3, 3) = negRtY;
	trans.block(3, 0, 3, 3) = zero;
	trans.block(3, 3, 3, 3) = Rt;
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
Matrix3d DynamicFormula::computeR_() {
	Matrix3d daOmega = toDaOmegaOrY(w);
	return R * daOmega;
}*/
Vector3d DynamicFormula::computey_() {
	return R * v;
}
VectorXd DynamicFormula::computelp() {
	VectorXd wv(6);
	wv.block(0, 0, 3, 1) = w;
	wv.block(3, 0, 3, 1) = v;
	return K * wv;
}
VectorXd DynamicFormula::computelp_() {
	Matrix3d Y = so3_ad(y);
	VectorXd tf = tsfs2tf(Y);
	Vector3d l = vec62Vec31(lp);
	Vector3d p = vec62Vec32(lp);
	Vector3d a = l.cross(w) + p.cross(v);
	Vector3d b = p.cross(w);
	VectorXd ab(6);
	ab.block(0, 0, 3, 1) = a;
	ab.block(3, 0, 3, 1) = b;
	return ab + tf;
}
Matrix3d DynamicFormula::s_GetRotaionMatrix(double angle, const Vector3d &axis)
{
	double x = axis(0);
	double y = axis(1);
	double z = axis(2);

#define  ZERO_TOL    1.0e-7
	double s = sin(angle);
	double c = cos(angle);

	Matrix3d matrix = Matrix3d::Identity();
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
	double mag = sqrt(x * x + y * y + z * z);
	mag = 1.0 / mag;
	x *= mag;
	y *= mag;
	z *= mag;
	double t = 1.0 - c;

	double tx = t * x;
	double ty = t * y;
	double tz = t * z;
	double sx = s * x;
	double sy = s * y;
	double sz = s * z;
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
Matrix3d DynamicFormula::computeNextR() {
	/*
	Matrix3d daomega = toDaOmegaOrY(w);
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
	Vector3d Rw = R * w;
	cout << "RW" << Rw << endl;
	double length = sqrt(Rw(0)*Rw(0) + Rw(1)*Rw(1) + Rw(2)*Rw(2));
	double angle = delta_t * length;
	Matrix3d delta_matrix = s_GetRotaionMatrix(angle, Rw);
	cout << "delta_matrix"<<delta_matrix << endl;
	return delta_matrix * R;*/
	/*
	cout << "q:" << q.coeffs() << endl;
	Vector3d Rw =R* w;
	cout << "RW" << Rw << endl;
	double length = sqrt(Rw(0)*Rw(0) + Rw(1)*Rw(1) + Rw(2)*Rw(2));
	theta = delta_t * length;
	double q0 = cos(delta_t*length / 2);
	double genhao = sin(delta_t*length / 2);
	double q1, q2, q3;
	if(length != 0){
		q1 = genhao * Rw(0) / length;
		q2 = genhao * Rw(1) / length;
		q3 = genhao * Rw(2) / length;
	}
	else {
		q1 = q2 = q3 = 0;
	}
	Quaternionf tempq(q0, q1, q2, q3);
	Matrix3d temR=tempq.toRotationMatrix();
	//q = tempq * q;//计算经过旋转后的orientation 
	//q.normalize();
	//return q.toRotationMatrix();
	return temR * R;*/
	/*
	cout << "q_old:" << q.coeffs() << endl;
	Vector3d Rw = R * w;
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
	Vector3d Rw = R *w;// 
	//double length = 1;sqrt(Rw(0)*Rw(0) + Rw(1)*Rw(1) + Rw(2)*Rw(2))
	Quaterniond q_w(0, Rw(0)* delta_t, Rw(1)*delta_t, Rw(2)* delta_t);
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
Vector3d DynamicFormula::computeNexty( Vector3d y_) {
	Vector3d temp_deltay = delta_t  * y_;
//	cout << "temp_deltay" << temp_deltay << endl;
//	cout << "(" << y(0) << "," << y(1) << "," << y(2) << ")" << endl;
	return y + temp_deltay;
}
VectorXd DynamicFormula::computeNextlp() {
	/*
	Vector3d l = lp.block(0, 0, 3, 1);
	Vector3d p = lp.block(3, 0, 3, 1);
	Vector3d l_ = lp_.block(0, 0, 3, 1);
	Vector3d p_ = lp_.block(3, 0, 3, 1);
	Quaternionf q_l(0, l_(0)*delta_t, l_(1) * delta_t, l_(2) * delta_t);
	q_l = q_l * l;
	q_l.w() = q_l.w()*0.5;
	q_l.x() = q_l.x()*0.5;
	q_l.y() = q_l.y()*0.5;
	q_l.z() = q_l.z()*0.5;
	cout << "q_l:" << q_l.coeffs() << endl;//计算经过旋转后的orientation 
	q.w() = q_l.w() + l.w();
	q.x() = q_l.x() + l.x();
	q.y() = q_l.y() + q.y();
	q.z() = q_l.z() + q.z();
	q.normalize();
	return q.toRotationMatrix();*/
	return lp +  lp_* delta_t ;
}
VectorXd DynamicFormula::computeNextwv() {
	MatrixXd Kinv = K.inverse();
	cout << "Kinv" << Kinv << endl;
	VectorXd res=Kinv * lp;
	//w = res.block(0, 0, 3, 1);
	//v = res.block(3, 0, 3, 1);
	return res;
}
Vector3d DynamicFormula::vec62Vec31(VectorXd wv) {//也可以用于wv
	return wv.block(0, 0, 3, 1);
}
Vector3d DynamicFormula::vec62Vec32(VectorXd wv) {
	return wv.block(3, 0, 3, 1);
}
Matrix3d DynamicFormula::so3_cay(Vector3d tempw) {
	Matrix3d daOmega = so3_ad(tempw);
	double wlength2 = tempw(0)*tempw(0) + tempw(1)*tempw(1) + tempw(2)*tempw(2);
	Matrix3d res = Matrix3d::Identity();
	double xishu = 4 / (4 + wlength2);
	res = res + xishu * (daOmega + daOmega * daOmega / 2);
	//cout << "so3cay" << res;
	return res;
}
Matrix4d DynamicFormula::se3_cay(VectorXd tempepsilon) {
	Vector3d tempw = tempepsilon.block(0, 0, 3, 1);//w
	Vector3d tempv = tempepsilon.block(3, 0, 3, 1);//v
	Matrix3d so3cay = so3_cay(delta_t*tempw);
	double wlength2 = tempw(0)*tempw(0) + tempw(1)*tempw(1) + tempw(2)*tempw(2);
	Matrix3d daOmega = so3_ad(tempw);
	double xishu = 2 / (4 + wlength2);
	Matrix3d tempmat = Matrix3d::Identity();
	tempmat = 2 * tempmat + daOmega;
	Matrix3d B = xishu * tempmat;
	Matrix4d res;
	res.block(0, 0, 3, 3) = so3cay;
	res.block(0, 3, 3, 1) = delta_t * B*tempv;
	//cout << "delta_t * B*v" << delta_t * B*tempv << endl;
	res(3, 0) = 0;
	res(3, 1) = 0;
	res(3, 2) = 0;
	res(3, 3) = 1;
	//cout << "se3_cay" << res << endl;
	return res;
}

MatrixXd DynamicFormula::se3_ad(VectorXd tempepsilon) {
	Vector3d w = tempepsilon.block(0, 0, 3, 1);//w
	Vector3d v = tempepsilon.block(3, 0, 3, 1);//v
	MatrixXd res(6, 6);
	res.block(0, 0, 3, 3) = so3_ad(w);
	res.block(0, 3, 3, 3) = Matrix3d::Zero();
	res.block(3, 0, 3, 3) = so3_ad(v);
	res.block(3, 3, 3, 3) = so3_ad(w);
	return res;
}
VectorXd DynamicFormula::se3_DEP(VectorXd epsilon_now, VectorXd epsilon_last, Matrix4d &gk) {
	Vector3d tempy = gk.block(0, 3, 3, 1);			
	Matrix3d Y = so3_ad(tempy);
	VectorXd tf = tsfs2tf(Y);
	MatrixXd ctln1 = se3_Ctln(delta_t * epsilon_now);
	MatrixXd ctln2 = se3_Ctln(-delta_t * epsilon_last);
	return ctln1*K*epsilon_now -
		ctln2*K*epsilon_last- delta_t* tf ;//
}

VectorXd DynamicFormula::Unconstr_Dyn(VectorXd epsilon_now, VectorXd epsilon_last, Matrix4d &gk) {
	VectorXd fepsilonk_est=se3_DEP(epsilon_now, epsilon_last, gk);//f epsilonk的估计值
	MatrixXd Jacobian(6, 6);
	Vector3d w = epsilon_now.block(0, 0, 3, 1);
	Vector3d v = epsilon_now.block(3, 0, 3, 1);
	Jacobian(0, 0) = K(0, 0) + delta_t * 0.5*w(2)*K(1, 0) - delta_t * 0.5*w(1)*K(2, 0);
	Jacobian(0, 1) = K(0, 1) + delta_t * 0.5*w(2)*K(1, 1) - delta_t * 0.5*w(1)*K(2, 1);
	Jacobian(0, 2) = K(0, 2) + delta_t * 0.5*w(2)*K(1, 2) - delta_t * 0.5*w(1)*K(2, 2);
	Jacobian(0, 3) = K(0, 3) + delta_t * 0.5*w(2)*K(1, 3) - delta_t * 0.5*w(1)*K(2, 3);
	Jacobian(0, 4) = K(0, 4) + delta_t * 0.5*w(2)*K(1, 4) - delta_t * 0.5*w(1)*K(2, 4);
	Jacobian(0, 5) = K(0, 5) + delta_t * 0.5*w(2)*K(1, 5) - delta_t * 0.5*w(1)*K(2, 5);

	Jacobian(1, 0) = -delta_t * 0.5*w(2)*K(0, 0) + K(1, 0) + delta_t * 0.5*w(0)*K(2, 0);
	Jacobian(1, 1) = -delta_t * 0.5*w(2)*K(0, 1) + K(1, 1) + delta_t * 0.5*w(0)*K(2, 1);
	Jacobian(1, 2) = -delta_t * 0.5*w(2)*K(0, 2) + K(1, 2) + delta_t * 0.5*w(0)*K(2, 2);
	Jacobian(1, 3) = -delta_t * 0.5*w(2)*K(0, 3) + K(1, 3) + delta_t * 0.5*w(0)*K(2, 3);
	Jacobian(1, 4) = -delta_t * 0.5*w(2)*K(0, 4) + K(1, 4) + delta_t * 0.5*w(0)*K(2, 4);
	Jacobian(1, 5) = -delta_t * 0.5*w(2)*K(0, 5) + K(1, 5) + delta_t * 0.5*w(0)*K(2, 5);

	Jacobian(2, 0) = delta_t * 0.5*w(1)*K(0, 0) - delta_t * 0.5*w(0)*K(1, 0) + K(2, 0);
	Jacobian(2, 1) = delta_t * 0.5*w(1)*K(0, 1) - delta_t * 0.5*w(0)*K(1, 1) + K(2, 1);
	Jacobian(2, 2) = delta_t * 0.5*w(1)*K(0, 2) - delta_t * 0.5*w(0)*K(1, 2) + K(2, 2);
	Jacobian(2, 3) = delta_t * 0.5*w(1)*K(0, 3) - delta_t * 0.5*w(0)*K(1, 3) + K(2, 3);
	Jacobian(2, 4) = delta_t * 0.5*w(1)*K(0, 4) - delta_t * 0.5*w(0)*K(1, 4) + K(2, 4);
	Jacobian(2, 5) = delta_t * 0.5*w(1)*K(0, 5) - delta_t * 0.5*w(0)*K(1, 5) + K(2, 5);

	Jacobian(3, 0) = delta_t * 0.5*v(2)*K(1, 0) - delta_t * 0.5*v(1)*K(2, 0) + K(3, 0) +
		delta_t * 0.5*w(2)*K(4, 0) - delta_t * 0.5*w(1)*K(5, 0);
	Jacobian(3, 1) = delta_t * 0.5*v(2)*K(1, 1) - delta_t * 0.5*v(1)*K(2, 1) + K(3, 1) +
		delta_t * 0.5*w(2)*K(4, 1) - delta_t * 0.5*w(1)*K(5, 1);
	Jacobian(3, 2) = delta_t * 0.5*v(2)*K(1, 2) - delta_t * 0.5*v(1)*K(2, 2) + K(3, 2) +
		delta_t * 0.5*w(2)*K(4, 2) - delta_t * 0.5*w(1)*K(5, 2);
	Jacobian(3, 3) = delta_t * 0.5*v(2)*K(1, 3) - delta_t * 0.5*v(1)*K(2, 3) + K(3, 3) +
		delta_t * 0.5*w(2)*K(4, 3) - delta_t * 0.5*w(1)*K(5, 3);
	Jacobian(3, 4) = delta_t * 0.5*v(2)*K(1, 4) - delta_t * 0.5*v(1)*K(2, 4) + K(3, 4) +
		delta_t * 0.5*w(2)*K(4, 4) - delta_t * 0.5*w(1)*K(5, 4);
	Jacobian(3, 5) = delta_t * 0.5*v(2)*K(1, 5) - delta_t * 0.5*v(1)*K(2, 5) + K(3, 5) +
		delta_t * 0.5*w(2)*K(4, 5) - delta_t * 0.5*w(1)*K(5, 5);

	Jacobian(4, 0) = -delta_t * 0.5*v(2)*K(0, 0) + delta_t * 0.5*v(0)*K(2, 0) -
		delta_t * 0.5*w(2)*K(3, 0) + K(4, 0) + delta_t * 0.5*w(0)*K(5, 0);
	Jacobian(4, 1) = -delta_t * 0.5*v(2)*K(0, 1) + delta_t * 0.5*v(0)*K(2, 1) -
		delta_t * 0.5*w(2)*K(3, 1) + K(4, 1) + delta_t * 0.5*w(0)*K(5, 1);
	Jacobian(4, 2) = -delta_t * 0.5*v(2)*K(0, 2) + delta_t * 0.5*v(0)*K(2, 2) -
		delta_t * 0.5*w(2)*K(3, 2) + K(4, 2) + delta_t * 0.5*w(0)*K(5, 2);
	Jacobian(4, 3) = -delta_t * 0.5*v(2)*K(0, 3) + delta_t * 0.5*v(0)*K(2, 3) -
		delta_t * 0.5*w(2)*K(3, 3) + K(4, 3) + delta_t * 0.5*w(0)*K(5, 3);
	Jacobian(4, 4) = -delta_t * 0.5*v(2)*K(0, 4) + delta_t * 0.5*v(0)*K(2, 4) -
		delta_t * 0.5*w(2)*K(3, 4) + K(4, 4) + delta_t * 0.5*w(0)*K(5, 4);
	Jacobian(4, 5) = -delta_t * 0.5*v(2)*K(0, 5) + delta_t * 0.5*v(0)*K(2, 5) -
		delta_t * 0.5*w(2)*K(3, 5) + K(4, 5) + delta_t * 0.5*w(0)*K(5, 5);

	Jacobian(5, 0) = delta_t * 0.5*v(1)*K(0, 0) - delta_t * 0.5*v(0)*K(1, 0) +
		delta_t * 0.5*w(1)*K(3, 0) - delta_t * 0.5*w(0)*K(4, 0) + K(5, 0);
	Jacobian(5, 1) = delta_t * 0.5*v(1)*K(0, 1) - delta_t * 0.5*v(0)*K(1, 1) +
		delta_t * 0.5*w(1)*K(3, 1) - delta_t * 0.5*w(0)*K(4, 1) + K(5, 1);
	Jacobian(5, 2) = delta_t * 0.5*v(1)*K(0, 2) - delta_t * 0.5*v(0)*K(1, 2) +
		delta_t * 0.5*w(1)*K(3, 2) - delta_t * 0.5*w(0)*K(4, 2) + K(5, 2);
	Jacobian(5, 3) = delta_t * 0.5*v(1)*K(0, 3) - delta_t * 0.5*v(0)*K(1, 3) +
		delta_t * 0.5*w(1)*K(3, 3) - delta_t * 0.5*w(0)*K(4, 3) + K(5, 3);
	Jacobian(5, 4) = delta_t * 0.5*v(1)*K(0, 4) - delta_t * 0.5*v(0)*K(1, 4) +
		delta_t * 0.5*w(1)*K(3, 4) - delta_t * 0.5*w(0)*K(4, 4) + K(5, 4);
	Jacobian(5, 5) = delta_t * 0.5*v(1)*K(0, 5) - delta_t * 0.5*v(0)*K(1, 5) +
		delta_t * 0.5*w(1)*K(3, 5) - delta_t * 0.5*w(0)*K(4, 5) + K(5, 5);
	VectorXd delta_epsilion=Jacobian.inverse() * fepsilonk_est;
	//cout << "delta_epsilion是是是=" << delta_epsilion << endl;
	return epsilon_now - delta_epsilion;
}
MatrixXd DynamicFormula::se3_Ctln(VectorXd tempepsilon) {
	MatrixXd identity(6, 6);
	identity.setIdentity();
	MatrixXd temp=0.5f*se3_ad(tempepsilon);
	return identity - temp;
}

void DynamicFormula::nextTime() {
	//lp_ = computelp_();
	//cout << "oldlp" << lp << endl;
	//cout << "nextlp_:" << lp_<< endl;
	//lp=computeNextlp();//lp
	//cout << "lp:" << lp << endl;
	//VectorXd tempwv= computeNextwv();
	//cout << "tempwv" << tempwv << endl;
	//cout << "w:"<<w(0)<<" " << w(1) << " " << w(2) << endl;
	//cout << "v:"<<v(0) << " " << v(1) << " " << v(2) << endl;
	Vector3d tempy = g.block(0, 3, 3, 1);
	Matrix3d Y = so3_ad(tempy);
	VectorXd tf = tsfs2tf(Y);
	VectorXd epsilon_last = epsilon;
	VectorXd delta_epsilon= delta_t * K.inverse()*(se3_ad(delta_t*epsilon_last)*K*epsilon_last +tf);
	VectorXd epsilon_now = epsilon_last + delta_epsilon;
	//cout << "epsilon_now 预估" << epsilon_now << endl;
	VectorXd res = se3_DEP(epsilon_now, epsilon_last, g);
	//cout << "偏差res=" << res << endl;
	//牛顿迭代法求解方程组
	int i = 0;
	double cancha =1e-15;
	do{//设置残差值
	
		epsilon_now=Unconstr_Dyn(epsilon_now, epsilon_last, g);
		//cout << "epsilon_now 迭代后的值" << epsilon_now << endl;
		res = se3_DEP(epsilon_now, epsilon_last, g);
		//cout << "偏差res=" << res << endl;
		i++;
	} while ((res(0) > cancha || res(1) > cancha || res(2) > cancha ||
		res(3) > cancha || res(4) > cancha || res(5) > cancha ||
		res(0) < -cancha || res(1) < -cancha || res(2) < -cancha ||
		res(3) < -cancha || res(4) < -cancha || res(5) < -cancha
		) && i < 50);
	//cout << "迭代了多少次？" << i << endl;
	Matrix4d se3cay= se3_cay(delta_t * epsilon_now);
	//cout << "se3cay==" << se3cay << endl;
	g= g * se3cay;
	//cout << "g" << g << endl;
	R = g.block(0, 0, 3, 3);
	y = g.block(0, 3, 3, 1);
	epsilon = epsilon_now;
	//Vector3d y_ = computey_();
	//cout << "y_" << y_ << endl;
	//y=computeNexty(y_);//y
	cout << "y!!!!!!!!!!!!!!:" << y << endl;
	//R=computeNextR();//R
	cout << "R:" << R << endl;
	//cout << "w:" << epsilon.block(0, 0, 3, 1) << endl;
	//cout << "v:" << epsilon.block(3, 0, 3, 1) << endl;
	//w = tempwv.block(0, 0, 3, 1);//w
	//v = tempwv.block(3, 0, 3, 1);//v
}
void DynamicFormula::set_tsfs(Vector3d ts,Vector3d fs) {

	//this->tsfs
}